
#include "Acts/EventData/SourceLink.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Definitions/Algebra.hpp"
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
#include "ActsExamples/Geometry/MaterialWiper.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/Io/Root/RootSimHitWriter.hpp"
#include "ActsExamples/Io/Root/RootSpacepointWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSeedWriter.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/Io/Root/RootTrackStatesWriter.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementWriter.hpp"
// Temporary
#include "ActsExamples/TruthTracking/TruthSeedingAlgorithm.hpp"

// Needed for detector construction
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"

#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
//  native units: mm, GeV, c=1, e=1

class MyDetectorElement : public Acts::DetectorElementBase
{
public:
  MyDetectorElement(std::shared_ptr<const Acts::Transform3> transform, std::shared_ptr<Acts::Surface> surface, double thickness)
      : Acts::DetectorElementBase(), mTransform(transform), mSurface(surface), mThickness(thickness)
  {
  }
  const Acts::Transform3 &transform(const Acts::GeometryContext &gctx) const { return *mTransform; }
  const Acts::Surface &surface() const { return *mSurface; }
  Acts::Surface &surface() { return *mSurface; }
  double thickness() const { return mThickness; }

private:
  std::shared_ptr<const Acts::Transform3> mTransform = nullptr;
  std::shared_ptr<Acts::Surface> mSurface = nullptr;
  double mThickness = 0.;
};

int main()
{
  // Logger
  Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;
  Acts::Logging::Level logLevelFatras = Acts::Logging::VERBOSE;
  Acts::Logging::Level logLevelDigi = Acts::Logging::VERBOSE;
  Acts::Logging::Level logLevelSeed = Acts::Logging::VERBOSE;
  Acts::Logging::Level logLevelFinder = Acts::Logging::VERBOSE;

  // Random number generator config
  auto rnd = std::make_shared<ActsExamples::RandomNumbers>(ActsExamples::RandomNumbers::Config({42}));

  // Particle gun generator config
  ActsExamples::ParametricParticleGenerator::Config genConfig;
  genConfig.etaUniform = true;
  genConfig.thetaMin = 0.30;
  genConfig.thetaMax = 0.40;
  genConfig.pMin = 3.0_GeV;
  genConfig.pMax = 3.1_GeV;
  ActsExamples::EventGenerator::Generator gen{
      std::make_shared<ActsExamples::FixedMultiplicityGenerator>(1),
      std::make_shared<ActsExamples::FixedVertexGenerator>(),
      std::make_shared<ActsExamples::ParametricParticleGenerator>(genConfig)};
  ActsExamples::EventGenerator::Config evgenConfig;
  evgenConfig.outputParticles = "particles";
  evgenConfig.generators = {gen};
  evgenConfig.randomNumbers = rnd;

  // Telescope
  // std::shared_ptr<const Acts::IMaterialDecorator> matDecorator = std::make_shared<const Acts::MaterialWiper>();
  // ActsExamples::Telescope::TelescopeDetector::Config telescopeConfig; // default telescope
  // ActsExamples::Telescope::TelescopeDetector telescope;
  // auto telescopePair = telescope.finalize(telescopeConfig, matDecorator); // std::pair<TrackingGeometryPtr, ContextDecorators>

  // Detector parameters
  std::vector<double> positions{{60, 90, 120, 150, 180, 210}}; // z-positions of sensitive layers
  double thickness{50_um};                                     // 
  double rMin{9_mm};
  double rMax{100_mm};

  // Create materials
  Acts::Material silicon = Acts::Material::fromMassDensity(9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);
  Acts::MaterialSlab matProp(silicon, thickness);
  const auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);

  // Construct the surfaces and layers
  std::vector<std::shared_ptr<MyDetectorElement>> detectorStore;
  Acts::LayerVector layVec;
  //const auto rBounds = std::make_shared<const Acts::RadialBounds>(rMin, rMax);  // <- for disk-like layers
  const auto pBounds = std::make_shared<const Acts::RectangleBounds>(rMax, rMax); // <- for square-like layers
  for (unsigned int i = 0; i < positions.size(); i++)
  {
    Acts::Translation3 trans(0, 0, positions[i]);
    Acts::Transform3 trafo(trans);
    // create surface
    // auto surface = Acts::Surface::makeShared<Acts::DiscSurface>(trafo, rMin, rMax); // <- for disk-like layers
    auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafo, pBounds); // <- for square-like layers
    surface->assignSurfaceMaterial(std::move(surfaceMaterial));
    auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(surface));
    // create layer
    // auto layer = Acts::DiscLayer::create(trafo, rBounds, std::move(surArray), 1._mm); // <- for disk-like layers
    auto layer = Acts::PlaneLayer::create(trafo, pBounds, std::move(surArray), 1._mm); // <- for square-like layers
    surface->associateLayer(*layer.get());
    layVec.push_back(layer);
    // create detector element
    auto detElement = std::make_shared<MyDetectorElement>(std::make_shared<const Acts::Transform3>(trafo), surface, 1._um);
    surface->assignDetectorElement(*detElement.get());
    // surface->assignGeometryId(i+256); // doesn't work?
    detectorStore.push_back(std::move(detElement));
  }

  // Create layer array
  Acts::GeometryContext context;
  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(lacConfig, Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(context, layVec, positions.front() - 2._mm, positions.back() + 2._mm, Acts::BinningType::arbitrary, Acts::BinningValue::binZ));

  // Build mother tracking volume
  Acts::Translation3 transVol(0, 0, (positions.front() + positions.back()) * 0.5);
  Acts::Transform3 trafoVol(transVol);
  auto length = positions.back() - positions.front();
  auto boundsVol = std::make_shared<const Acts::CuboidVolumeBounds>(rMax + 5._mm, rMax + 5._mm, length + 10._mm);
  // auto boundsVol = std::make_shared<const Acts::CylinderVolumeBounds>(rMin - 5._mm, rMax + 5._mm, length + 10._mm);
  auto trackVolume = Acts::TrackingVolume::create(trafoVol, boundsVol, nullptr, std::move(layArr), nullptr, {}, "Telescope");

  // Build tracking geometry
  auto trackingGeometry = std::make_shared<Acts::TrackingGeometry>(trackVolume);

  // Check geometry
  // Why duplicated layers?
  // Why geometryIds are not set?
  const Acts::TrackingVolume *highestTrackingVolume = trackingGeometry->highestTrackingVolume();
  // const Acts::TrackingVolume* highestTrackingVolume = telescopePair.first->highestTrackingVolume();
  printf("volumeId = %d\n", highestTrackingVolume->geometryId().value());
  const Acts::LayerArray *confinedLayers = highestTrackingVolume->confinedLayers();
  for (const auto &layer : confinedLayers->arrayObjects())
  {
    printf("  layerId = %d, thickness = %f\n", layer->geometryId().value(), layer->thickness());
    if (!layer->surfaceArray())
      continue;
    for (const auto &surface : layer->surfaceArray()->surfaces())
    {
      printf("    surfaceId = %d\n", surface->geometryId().value());
    }
  }

  // Fatras config
  ActsExamples::FatrasSimulation::Config fatrasConfig;
  fatrasConfig.inputParticles = evgenConfig.outputParticles;
  fatrasConfig.outputParticlesInitial = "initial";
  fatrasConfig.outputParticlesFinal = "final";
  fatrasConfig.outputSimHits = "simhits";
  // fatrasConfig.trackingGeometry = telescopePair.first;
  fatrasConfig.trackingGeometry = trackingGeometry;
  fatrasConfig.magneticField = std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, 1._T));
  fatrasConfig.randomNumbers = rnd;

  // Digitization config
  bool doMerge = 0;
  double mergeSigma = 0.001; // mm
  bool mergeCommonCorner = 0;
  ActsExamples::DigitizationConfig digiCfg(doMerge, mergeSigma, mergeCommonCorner);
  digiCfg.inputSimHits = fatrasConfig.outputSimHits;
  digiCfg.randomNumbers = rnd;
  digiCfg.trackingGeometry = trackingGeometry;
  digiCfg.digitizationConfigs = ActsExamples::readDigiConfigFromJson("digi-smearing-config.json");

  // Create space points
  ActsExamples::SpacePointMaker::Config spConfig;
  spConfig.inputSourceLinks = digiCfg.outputSourceLinks;
  spConfig.inputMeasurements = digiCfg.outputMeasurements;
  spConfig.trackingGeometry = trackingGeometry;
  spConfig.outputSpacePoints = "spacepoints";
  spConfig.geometrySelection = {Acts::GeometryIdentifier{}};

  // ActsExamples::TruthSeedingAlgorithm::Config truthSeedingCfg;
  // truthSeedingCfg.inputParticles = evgenConfig.outputParticles;
  // truthSeedingCfg.inputSpacePoints = {spConfig.outputSpacePoints};
  // truthSeedingCfg.inputMeasurementParticlesMap = digiCfg.outputMeasurementParticlesMap;
  // truthSeedingCfg.outputSeeds = "seeds";
  // truthSeedingCfg.outputParticles = "seeded_particles";
  // truthSeedingCfg.outputProtoTracks = "prototracks";
  // truthSeedingCfg.deltaRMin = 1._mm;
  // truthSeedingCfg.deltaRMax = 100._mm;

  // Seeding
  ActsExamples::SeedingAlgorithm::Config seedingCfg;
  seedingCfg.inputSpacePoints = {spConfig.outputSpacePoints};
  seedingCfg.outputSeeds = "seeds";
  seedingCfg.gridOptions.bFieldInZ = 1_T;
  seedingCfg.gridConfig.minPt = 0.5_GeV;
  seedingCfg.gridConfig.zMin = positions[0] - 1_mm;
  seedingCfg.gridConfig.zMax = positions[2] + 1_mm;
  seedingCfg.gridConfig.rMax = rMax;
  seedingCfg.gridConfig.deltaRMax = rMax - rMin;
  seedingCfg.gridConfig.cotThetaMax = 7.0;
  seedingCfg.gridConfig.impactMax = 8_mm;
  seedingCfg.seedFilterConfig.maxSeedsPerSpM = 1;
  seedingCfg.seedFinderOptions.bFieldInZ        = seedingCfg.gridOptions.bFieldInZ;
  seedingCfg.seedFinderConfig.impactMax         = seedingCfg.gridConfig.impactMax;
  seedingCfg.seedFinderConfig.rMax              = seedingCfg.gridConfig.rMax;
  seedingCfg.seedFinderConfig.deltaRMax         = seedingCfg.gridConfig.deltaRMax;
  seedingCfg.seedFinderConfig.deltaRMaxTopSP    = seedingCfg.gridConfig.deltaRMax;
  seedingCfg.seedFinderConfig.deltaRMaxBottomSP = seedingCfg.gridConfig.deltaRMax;
  seedingCfg.seedFinderConfig.zMin              = seedingCfg.gridConfig.zMin;
  seedingCfg.seedFinderConfig.zMax              = seedingCfg.gridConfig.zMax;
  seedingCfg.seedFinderConfig.cotThetaMax       = seedingCfg.gridConfig.cotThetaMax;
  seedingCfg.seedFinderConfig.minPt             = seedingCfg.gridConfig.minPt;
  seedingCfg.seedFinderConfig.maxSeedsPerSpM    = seedingCfg.seedFilterConfig.maxSeedsPerSpM;
  seedingCfg.seedFinderConfig.deltaZMax = positions[2] - positions[0];
  seedingCfg.seedFinderConfig.rMinMiddle = rMin;
  seedingCfg.seedFinderConfig.rMaxMiddle = rMax;

  // Parameter Estimation from seeds
  ActsExamples::TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
  paramsEstimationCfg.inputSeeds = seedingCfg.outputSeeds;
  paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
  paramsEstimationCfg.trackingGeometry = trackingGeometry;
  paramsEstimationCfg.magneticField = fatrasConfig.magneticField;

  // Track finding
  ActsExamples::TrackFindingAlgorithm::Config trackFindingCfg;
  trackFindingCfg.inputMeasurements = digiCfg.outputMeasurements;
  trackFindingCfg.inputSourceLinks = digiCfg.outputSourceLinks;
  trackFindingCfg.inputInitialTrackParameters = paramsEstimationCfg.outputTrackParameters;
  trackFindingCfg.outputTracks = "tracks";
  trackFindingCfg.measurementSelectorCfg = {{Acts::GeometryIdentifier(),{{}, {std::numeric_limits<double>::max()}, {1u}}}}; //chi2cut, numberOfMeasurementsPerSurface
  trackFindingCfg.findTracks = ActsExamples::TrackFindingAlgorithm::makeTrackFinderFunction(
    trackingGeometry, fatrasConfig.magneticField,*Acts::getDefaultLogger("TrackFinder", logLevelFinder));

  // write track states from CKF
  ActsExamples::RootTrackStatesWriter::Config trackStatesWriter;
  trackStatesWriter.inputTracks = trackFindingCfg.outputTracks;
  trackStatesWriter.inputParticles = evgenConfig.outputParticles;
  trackStatesWriter.inputSimHits = fatrasConfig.outputSimHits;
  trackStatesWriter.inputMeasurementParticlesMap = digiCfg.outputMeasurementParticlesMap;
  trackStatesWriter.inputMeasurementSimHitsMap = digiCfg.outputMeasurementSimHitsMap;
  trackStatesWriter.filePath = "trackstates_ckf.root";
  trackStatesWriter.treeName = "trackstates";

  // SimhitWriter config
  ActsExamples::RootSimHitWriter::Config simhitWriterConfig;
  simhitWriterConfig.inputSimHits = fatrasConfig.outputSimHits;
  simhitWriterConfig.filePath = "hits.root";

  // Measurement writer
  ActsExamples::CsvMeasurementWriter::Config measWriterCfg;
  measWriterCfg.inputMeasurements = digiCfg.outputMeasurements;
  measWriterCfg.inputMeasurementSimHitsMap = digiCfg.outputMeasurementSimHitsMap;
  measWriterCfg.outputDir = "measurements";

  // SpacepointWriter config
  ActsExamples::RootSpacepointWriter::Config spWriterConfig;
  spWriterConfig.inputSpacepoints = spConfig.outputSpacePoints;
  spWriterConfig.filePath = "spacepoints.root";

  // Seed writer config
  ActsExamples::CsvSeedWriter::Config seedWriterConfig;
  seedWriterConfig.inputSimHits = fatrasConfig.outputSimHits;
  seedWriterConfig.inputSimSeeds = seedingCfg.outputSeeds;
  seedWriterConfig.inputTrackParameters = paramsEstimationCfg.outputTrackParameters;
  seedWriterConfig.inputMeasurementParticlesMap = digiCfg.outputMeasurementParticlesMap;
  seedWriterConfig.inputMeasurementSimHitsMap = digiCfg.outputMeasurementSimHitsMap;
  seedWriterConfig.outputDir = "seeds";

  // Sequencer config
  ActsExamples::Sequencer::Config sequencerConfig;
  sequencerConfig.numThreads = 1;
  sequencerConfig.events = 10;
  sequencerConfig.logLevel = logLevel;

  // Start sequencer
  ActsExamples::Sequencer sequencer(sequencerConfig);
  sequencer.addReader(std::make_shared<ActsExamples::EventGenerator>(evgenConfig, logLevel));
  sequencer.addElement(std::make_shared<ActsExamples::FatrasSimulation>(fatrasConfig, logLevelFatras));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::DigitizationAlgorithm>(digiCfg, logLevelDigi));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::SpacePointMaker>(spConfig, logLevel));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::SeedingAlgorithm>(seedingCfg, logLevelSeed));
  //sequencer.addAlgorithm(std::make_shared<ActsExamples::TruthSeedingAlgorithm>(truthSeedingCfg, logLevelSeed));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::TrackParamsEstimationAlgorithm>(paramsEstimationCfg, logLevel));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::TrackFindingAlgorithm>(trackFindingCfg, logLevelFinder));
  sequencer.addWriter(std::make_shared<ActsExamples::RootSimHitWriter>(simhitWriterConfig, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::CsvMeasurementWriter>(measWriterCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootSpacepointWriter>(spWriterConfig, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::CsvSeedWriter>(seedWriterConfig, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootTrackStatesWriter>(trackStatesWriter, logLevel));

  sequencer.run();

  return 0;
}
