#ifndef tracker
#define tracker
#include "tracker_config.h"

class MyDetectorElement : public Acts::DetectorElementBase {
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

std::vector<std::shared_ptr<MyDetectorElement>> detectorStore;
Acts::GeometryContext gctx;

Acts::TrackingGeometry* CreateTrackingGeometry(){
  // Create materials
  double radLenSilicon = 9.370_cm;
  double radLenAluminium = 8.897_cm;
  Acts::Material silicon   = Acts::Material::fromMassDensity(9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);
  Acts::Material aluminium = Acts::Material::fromMassDensity(8.897_cm, 39.70_cm, 26.9815, 13, 2.699_g / 1_cm3);

  Acts::MaterialSlab matProp(silicon, thickness);
  const auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);

  // Construct the surfaces and layers
  Acts::LayerVector layVec;

  if (1) {
    double eps = 1e-10;
    double radLenFractionROC = 0.25;
    double thicknessROC = radLenFractionROC*radLenAluminium;
    double zDrift = 163_cm;
    double zROC = zDrift + thicknessROC/2.;
    double rMinROC = 35_cm;
    double rMaxROC = 140_cm;

    double radLenFractionFrame = 0.75;
    double thicknessFrame = radLenFractionFrame*radLenAluminium;

    double zFrame = zDrift + thicknessROC + thicknessFrame/2. + 1._mm;
    double zFrameCircum1 = zFrame;
    double rMinFrameCircum1 = 35_cm;
    double rMaxFrameCircum1 = 42_cm;

    double zFrameCircum2 = zFrame;
    double rMinFrameCircum2 = 120_cm;
    double rMaxFrameCircum2 = 140_cm;

    double zFrameRadial = zFrame;
    double halfYFrameRadial = 6.3_cm;
    double halfXFrameRadial = (sqrt(rMinFrameCircum2*rMinFrameCircum2-halfYFrameRadial*halfYFrameRadial)-rMaxFrameCircum1)/2.-2*eps;
    double cXFrameRadial0 = rMaxFrameCircum1 + halfXFrameRadial + eps;
    const int nSectors = 12;
    double phi[nSectors];
    for (int i=0;i<nSectors;i++){
      phi[i] = (15.+30.*i)/180.*M_PI;
    }
    
    Acts::MaterialSlab matPropROC(aluminium, thicknessROC);
    const auto surfaceMatROC = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matPropROC);
    Acts::Translation3 trans(0, 0, zROC);
    Acts::Transform3 trafo(trans);
    const auto rBoundsROC = std::make_shared<const Acts::RadialBounds>(rMinROC, rMaxROC);
    auto surface = Acts::Surface::makeShared<Acts::DiscSurface>(trafo, rMinROC, rMaxROC);
    surface->assignSurfaceMaterial(std::move(surfaceMatROC));
    auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(surface));
    auto layer = Acts::DiscLayer::create(trafo, rBoundsROC, std::move(surArray), thicknessROC, nullptr, Acts::active);
    surface->associateLayer(*layer.get());
    layVec.push_back(layer);

    if (0) { // create detector element to check geometry
    auto detElement = std::make_shared<MyDetectorElement>(std::make_shared<const Acts::Transform3>(trafo), surface, thicknessROC);
    surface->assignDetectorElement(*detElement.get());
    detectorStore.push_back(std::move(detElement));
    }
    
    if (1) {
    Acts::MaterialSlab matPropFrame(aluminium, thicknessFrame);
    const auto surfaceMatFrame = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matPropFrame);
    std::vector<std::shared_ptr<const Acts::Surface>> vSurfaceFrame;
    std::vector<Acts::Transform3> vTrafoFrame;

    Acts::Translation3 transFrameCircum1(0, 0, zFrameCircum1);
    Acts::Translation3 transFrameCircum2(0, 0, zFrameCircum2);
    Acts::Transform3 trafoFrameCircum1(transFrameCircum1);
    Acts::Transform3 trafoFrameCircum2(transFrameCircum2);

    const auto pBoundsFrameRadial = std::make_shared<const Acts::RectangleBounds>(halfXFrameRadial, halfYFrameRadial);
    for (int iSector = 0; iSector < nSectors; iSector++) {
      Acts::Transform3 trafoFrameRadial;
      trafoFrameRadial.setIdentity();
      trafoFrameRadial.rotate(Eigen::AngleAxisd(phi[iSector], Acts::Vector3(0, 0, 1)));
      trafoFrameRadial.translate(Acts::Vector3(cXFrameRadial0, 0, zFrameRadial));
      vTrafoFrame.push_back(trafoFrameRadial);
      auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafoFrameRadial, pBoundsFrameRadial);
      surface->assignSurfaceMaterial(std::move(surfaceMatFrame));
      vSurfaceFrame.push_back(surface);

      auto bounds1 = std::make_shared<const Acts::RadialBounds>(rMinFrameCircum1, rMaxFrameCircum1, M_PI/12.-eps, Acts::detail::radian_sym(iSector*M_PI/6));
      auto bounds2 = std::make_shared<const Acts::RadialBounds>(rMinFrameCircum2, rMaxFrameCircum2, M_PI/12.-eps, Acts::detail::radian_sym(iSector*M_PI/6));
      auto surfaceFrameCircum1 = Acts::Surface::makeShared<Acts::DiscSurface>(trafoFrameCircum1, bounds1);
      auto surfaceFrameCircum2 = Acts::Surface::makeShared<Acts::DiscSurface>(trafoFrameCircum2, bounds2);
      surfaceFrameCircum1->assignSurfaceMaterial(std::move(surfaceMatFrame));
      surfaceFrameCircum2->assignSurfaceMaterial(std::move(surfaceMatFrame));
      vSurfaceFrame.push_back(surfaceFrameCircum1);
      vSurfaceFrame.push_back(surfaceFrameCircum2);
      vTrafoFrame.push_back(trafoFrameCircum1);
      vTrafoFrame.push_back(trafoFrameCircum2);
    }

    Acts::SurfaceArrayCreator::Config sacConfig;
    Acts::SurfaceArrayCreator surfaceArrayCreator(sacConfig, Acts::getDefaultLogger("SurfaceArrayCreator", Acts::Logging::VERBOSE));
    const auto rBoundsFrame = std::make_shared<const Acts::RadialBounds>(rMinFrameCircum1, rMaxFrameCircum2);
    Acts::Translation3 transFrame(0., 0., zFrameRadial);
    Acts::Transform3 trafoFrame(transFrame);
    auto layerFrame = Acts::DiscLayer::create(trafoFrame, rBoundsFrame, std::move(surfaceArrayCreator.surfaceArrayOnDisc(gctx, vSurfaceFrame, 3, 12)), thicknessFrame, nullptr, Acts::active);

    int i=0;
    for (auto& surface : vSurfaceFrame) {
      auto mutableSurface = const_cast<Acts::Surface*>(surface.get());
      mutableSurface->associateLayer(*layerFrame.get());
      if (0) { // create detector element to check geometry
      auto detElement = std::make_shared<MyDetectorElement>(std::make_shared<const Acts::Transform3>(vTrafoFrame[i]), std::const_pointer_cast<Acts::Surface>(surface), thicknessFrame);
      mutableSurface->assignDetectorElement(*detElement.get());
      detectorStore.push_back(std::move(detElement));
      }
      i++;
    }
    layVec.push_back(layerFrame);
    }
  }


  // const auto rBounds = std::make_shared<const Acts::RadialBounds>(rMin, rMax);  // <- for disk-like layers
  const auto pBounds = std::make_shared<const Acts::RectangleBounds>(rMax, rMax); // <- for square-like layers
  for (unsigned int i = 0; i < positions.size(); i++) {
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
    auto detElement = std::make_shared<MyDetectorElement>(std::make_shared<const Acts::Transform3>(trafo), surface, thickness);
    surface->assignDetectorElement(*detElement.get());
    detectorStore.push_back(std::move(detElement));
  }

  // Create layer array
  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(lacConfig, Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::VERBOSE));
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(gctx, layVec, 1500, positions.back() + 2._mm, Acts::BinningType::arbitrary, Acts::BinningValue::binZ));

  // Build mother tracking volume
  Acts::Translation3 transVol(0, 0, 0);
  Acts::Transform3 trafoVol(transVol);
  auto boundsVol = std::make_shared<const Acts::CuboidVolumeBounds>(2000._mm, 2000._mm, 3100._mm); // <- for square-like layers
  auto trackVolume = std::make_shared<Acts::TrackingVolume>(trafoVol, boundsVol, nullptr, std::move(layArr), nullptr, Acts::MutableTrackingVolumeVector{}, "Telescope");
  // Build tracking geometry
  auto trackingGeometry = new Acts::TrackingGeometry(trackVolume);

  // Check geometry
  // Why geometryIds are not set?
  const Acts::TrackingVolume *highestTrackingVolume = trackingGeometry->highestTrackingVolume();
  printf("volumeId = %d\n", highestTrackingVolume->geometryId().value());
  const Acts::LayerArray *confinedLayers = highestTrackingVolume->confinedLayers();
  for (const auto &layer : confinedLayers->arrayObjects())  {
    printf("  layerId = %d, thickness = %f type = %d\n", layer->geometryId().value(), layer->thickness(), layer->layerType());
    if (layer->layerType()==-1) {
      const Acts::NavigationLayer* nlayer = dynamic_cast<const Acts::NavigationLayer*>(layer.get());
      printf(" navigation %p\n", nlayer);
    } else if (layer->layerType()==0) {
      const Acts::DiscLayer* player = dynamic_cast<const Acts::DiscLayer*>(layer.get());
      printf(" disk %p\n", player);
    } else {
      const Acts::PlaneLayer* player = dynamic_cast<const Acts::PlaneLayer*>(layer.get());
      printf(" plane %p\n", player);
    }
    if (!layer->surfaceArray())
      continue;
    for (const auto &surface : layer->surfaceArray()->surfaces())
    {
      printf("    surfaceId = %d\n", surface->geometryId().value());
    }
  }

  return trackingGeometry;
}

#endif