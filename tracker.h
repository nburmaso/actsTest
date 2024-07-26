#ifndef tracker
#define tracker
#include "tracker_config.h"
#include "Acts/Definitions/Units.hpp"

using Acts::UnitConstants::cm;

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

Acts::TrackingGeometry* CreateTrackingGeometry(bool addROC = 1, bool addFlange = 1){
  // Create materials
  Acts::Material silicon   = Acts::Material::fromMassDensity(9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);
  Acts::Material aluminium = Acts::Material::fromMassDensity(8.897_cm, 39.70_cm, 26.9815, 13, 2.699_g / 1_cm3);

  Acts::MaterialSlab matProp(silicon, thickness*cm);
  const auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);

  // Construct the surfaces and layers
  Acts::LayerVector layVec;
 

  if (addROC) {
    double phi[nSectors];
    for (int i=0;i<nSectors;i++){
      phi[i] = (15.+30.*i)/180.*M_PI;
    }
    
    Acts::MaterialSlab matPropROC(aluminium, thicknessROC*cm);
    const auto surfaceMatROC = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matPropROC);
    Acts::Translation3 trans(0, 0, zROC*cm);
    Acts::Transform3 trafo(trans);
    const auto rBoundsROC = std::make_shared<const Acts::RadialBounds>(rMinROC*cm, rMaxROC*cm);
    auto surface = Acts::Surface::makeShared<Acts::DiscSurface>(trafo, rMinROC*cm, rMaxROC*cm);
    surface->assignSurfaceMaterial(std::move(surfaceMatROC));
    auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(surface));
    auto layer = Acts::DiscLayer::create(trafo, rBoundsROC, std::move(surArray), thicknessROC*cm, nullptr, Acts::active);
    surface->associateLayer(*layer.get());
    layVec.push_back(layer);

    if (addFlange) {
      Acts::MaterialSlab matPropFrame(aluminium, thicknessFrame*cm);
      const auto surfaceMatFrame = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matPropFrame);
      std::vector<std::shared_ptr<const Acts::Surface>> vSurfaceFrame;

      Acts::Translation3 transFrameCircum1(0, 0, zFrameCircum1*cm);
      Acts::Translation3 transFrameCircum2(0, 0, zFrameCircum2*cm);
      Acts::Transform3 trafoFrameCircum1(transFrameCircum1);
      Acts::Transform3 trafoFrameCircum2(transFrameCircum2);

      const auto pBoundsFrameRadial = std::make_shared<const Acts::RectangleBounds>(halfXFrameRadial*cm, halfYFrameRadial*cm);
      for (int iSector = 0; iSector < nSectors; iSector++) {
        Acts::Transform3 trafoFrameRadial;
        trafoFrameRadial.setIdentity();
        trafoFrameRadial.rotate(Eigen::AngleAxisd(phi[iSector], Acts::Vector3(0, 0, 1)));
        trafoFrameRadial.translate(Acts::Vector3(cXFrameRadial0*cm, 0, zFrameRadial*cm));
        auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafoFrameRadial, pBoundsFrameRadial);
        surface->assignSurfaceMaterial(std::move(surfaceMatFrame));
        vSurfaceFrame.push_back(surface);

        auto bounds1 = std::make_shared<const Acts::RadialBounds>(rMinFrameCircum1*cm, rMaxFrameCircum1*cm, M_PI/12.-eps, Acts::detail::radian_sym(iSector*M_PI/6));
        auto bounds2 = std::make_shared<const Acts::RadialBounds>(rMinFrameCircum2*cm, rMaxFrameCircum2*cm, M_PI/12.-eps, Acts::detail::radian_sym(iSector*M_PI/6));
        auto surfaceFrameCircum1 = Acts::Surface::makeShared<Acts::DiscSurface>(trafoFrameCircum1, bounds1);
        auto surfaceFrameCircum2 = Acts::Surface::makeShared<Acts::DiscSurface>(trafoFrameCircum2, bounds2);
        surfaceFrameCircum1->assignSurfaceMaterial(std::move(surfaceMatFrame));
        surfaceFrameCircum2->assignSurfaceMaterial(std::move(surfaceMatFrame));
        vSurfaceFrame.push_back(surfaceFrameCircum1);
        vSurfaceFrame.push_back(surfaceFrameCircum2);
      }

      Acts::SurfaceArrayCreator::Config sacConfig;
      Acts::SurfaceArrayCreator surfaceArrayCreator(sacConfig, Acts::getDefaultLogger("SurfaceArrayCreator", Acts::Logging::INFO));
      const auto rBoundsFrame = std::make_shared<const Acts::RadialBounds>(rMinFrameCircum1*cm, rMaxFrameCircum2*cm);
      Acts::Translation3 transFrame(0., 0., zFrameRadial*cm);
      Acts::Transform3 trafoFrame(transFrame);
      auto layerFrame = Acts::DiscLayer::create(trafoFrame, rBoundsFrame, std::move(surfaceArrayCreator.surfaceArrayOnDisc(gctx, vSurfaceFrame, 3, 12)), thicknessFrame*cm, nullptr, Acts::active);

      for (auto& surface : vSurfaceFrame) {
        auto mutableSurface = const_cast<Acts::Surface*>(surface.get());
        mutableSurface->associateLayer(*layerFrame.get());
      }
      layVec.push_back(layerFrame);
    }
  }


  // const auto rBounds = std::make_shared<const Acts::RadialBounds>(rMinStation, rMaxStation);  // <- for disk-like layers
  const auto pBounds = std::make_shared<const Acts::RectangleBounds>(rMaxStation*cm, rMaxStation*cm); // <- for square-like layers
  for (unsigned int i = 0; i < positions.size(); i++) {
    Acts::Translation3 trans(0, 0, positions[i]*cm);
    Acts::Transform3 trafo(trans);
    // create surface
    // auto surface = Acts::Surface::makeShared<Acts::DiscSurface>(trafo, rMinStation, rMaxStation); // <- for disk-like layers
    auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafo, pBounds); // <- for square-like layers
    surface->assignSurfaceMaterial(std::move(surfaceMaterial));
    auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(surface));
    // create layer
    // auto layer = Acts::DiscLayer::create(trafo, rBounds, std::move(surArray), 1._mm); // <- for disk-like layers
    auto layer = Acts::PlaneLayer::create(trafo, pBounds, std::move(surArray), 1._mm); // <- for square-like layers
    surface->associateLayer(*layer.get());
    layVec.push_back(layer);
    // create detector element
    auto detElement = std::make_shared<MyDetectorElement>(std::make_shared<const Acts::Transform3>(trafo), surface, thickness*cm);
    surface->assignDetectorElement(*detElement.get());
    detectorStore.push_back(std::move(detElement));
  }

  // Create layer array
  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(lacConfig, Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(gctx, layVec, 1500_mm, positions.back()*cm + 2._mm, Acts::BinningType::arbitrary, Acts::BinningValue::binZ));

  // Build mother tracking volume
  Acts::Translation3 transVol(0, 0, 0);
  Acts::Transform3 trafoVol(transVol);
  auto boundsVol = std::make_shared<const Acts::CuboidVolumeBounds>(2000._mm, 2000._mm, 3100._mm); // <- for square-like layers
  auto trackVolume = std::make_shared<Acts::TrackingVolume>(trafoVol, boundsVol, nullptr, std::move(layArr), nullptr, Acts::MutableTrackingVolumeVector{}, "Telescope");
  // Build tracking geometry
  auto trackingGeometry = new Acts::TrackingGeometry(trackVolume);

  // Check geometry
  // Why geometryIds are not set?
  bool print_surface_info = 0;
  if (print_surface_info) {
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
    } //for layers
  } // print_surface_info
  return trackingGeometry;
}

#endif