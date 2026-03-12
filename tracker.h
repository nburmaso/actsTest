#ifndef tracker
#define tracker
#include "tracker_config.h"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/NavigationLayer.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Utilities/Helpers.hpp"
//#include "Acts/Plugins/Geant4/Geant4Converters.hpp"

//#include "Geant4/G4Material.hh"
//#include "Geant4/G4NistManager.hh"


using Acts::UnitConstants::cm;
using namespace Acts::UnitLiterals;

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

Acts::TrackingGeometry* CreateTrackingGeometry(bool addROC = 0, bool addFlange = 0){
  auto gctx = Acts::GeometryContext();

  auto aluminium = Acts::Material::fromMassDensity(8.897_cm, 39.70_cm, 26.9815, 13, 2.699_g / 1_cm3);
  auto silicon = Acts::Material::fromMassDensity(9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);
  auto vacuum = Acts::Material::fromMassDensity(1e+10_cm, 1e+10_cm, 1.0000, 1, 1e-10_g / 1_cm3);
  auto volumeMaterial = std::make_shared<Acts::HomogeneousVolumeMaterial>(vacuum);

  std::vector<float> layerBinEdges;
  std::vector<std::pair<std::shared_ptr<const Acts::Layer>, Acts::Vector3>> layerOrderVec;

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
    layerOrderVec.emplace_back(layer, layer->referencePosition(gctx, Acts::AxisDirection::AxisZ));
    layerBinEdges.push_back(zROC*cm - thicknessROC/2.*cm);
    layerBinEdges.push_back(zROC*cm + thicknessROC/2.*cm);

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
      layerOrderVec.emplace_back(layerFrame, layerFrame->referencePosition(gctx, Acts::AxisDirection::AxisZ));
      layerBinEdges.push_back(zFrameRadial*cm + thicknessFrame/2.*cm);
    }
  }

  // add first navigation layer
  Acts::Transform3 tr1 = Acts::Transform3::Identity();
  tr1.translate(Acts::Vector3(0, 0, positions[0]*cm - 20));
  const auto pBounds = std::make_shared<const Acts::RectangleBounds>(rMaxStation*cm, rMaxStation*cm);
  auto layer = Acts::PlaneLayer::create(tr1, pBounds, nullptr, 1_mm, nullptr, Acts::navigation);
  layerOrderVec.emplace_back(layer, layer->referencePosition(gctx, Acts::AxisDirection::AxisZ));
  layerBinEdges.push_back(positions[0]*cm - 30); // TODO adjust bin edges
  layerBinEdges.push_back(positions[0]*cm - 10); // TODO adjust bin edges

  for (int i = 0; i < positions.size(); ++i) {
    double rmin = layerRMin[i] * cm;
    double rmax = layerRMax[i] * cm;
    layerBinEdges.push_back(i < positions.size()-1 ? (positions[i] + positions[i+1])/2.*cm : positions.back()*cm + 1);
    double thickness = layerThickness[i] * cm;
    // create surface material
    auto matProp = Acts::MaterialSlab(silicon, thickness);
    auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);

    if (layerType[i]==2) {
      const auto pBounds = std::make_shared<const Acts::RectangleBounds>(rmax, rmax);
      Acts::Transform3 trafo = Acts::Transform3::Identity();
      trafo.translate(Acts::Vector3(0, 0, positions[i]*cm));
      auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafo, pBounds);
      surface->assignSurfaceMaterial(std::move(surfaceMaterial));
      auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(surface));
      auto layer = Acts::PlaneLayer::create(trafo, pBounds, std::move(surArray), thickness, nullptr, Acts::active);
      surface->associateLayer(*layer.get());
      layerOrderVec.emplace_back(layer, layer->referencePosition(gctx, Acts::AxisDirection::AxisZ));
      auto detElement = std::make_shared<MyDetectorElement>(std::make_shared<const Acts::Transform3>(trafo), surface, thickness);
      surface->assignDetectorElement(*detElement.get());
      detectorStore.push_back(std::move(detElement));
      continue;
    }

    double angleRot = 0;
    if (layerType[i]==4) angleRot = 0;
    if (layerType[i]==5) angleRot = angleRotDeg*M_PI/180.;
    if (layerType[i]==6) angleRot =-angleRotDeg*M_PI/180.;
    double rc = 0.5 * (rmax + rmin);   
    double hl = 0.5 * (rmax - rmin) / cos(angleRot);
    double hw = layerStripWidth/2. * cm;
    const auto sBounds = std::make_shared<const Acts::RectangleBounds>(hw, hl);
    //const auto sBounds = std::make_shared<const Acts::LineBounds>(hw, hl); // for strawSurface
    int nTubes = numberOfTubes[i];
    std::vector<std::shared_ptr<const Acts::Surface>> vSurfaces;
    for (int iTube = 0; iTube < nTubes; iTube++) {
      double phi = layerAngle[i] + 2*M_PI/nTubes*(iTube+0.5);
      double cosp = cos(phi);
      double sinp = sin(phi);
      auto stransform = Acts::Transform3::Identity();
      stransform.translate(Acts::Vector3(rc*cosp, rc*sinp, positions[i] * cm));
      stransform.rotate(Eigen::AngleAxisd(M_PI/2 + phi + angleRot, Acts::Vector3(0, 0, 1)));
      // stransform.rotate(Eigen::AngleAxisd(M_PI/2, Acts::Vector3(1, 0, 0))); // for strawSurface
      auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(stransform, sBounds);
      //auto surface = Acts::Surface::makeShared<Acts::StrawSurface>(stransform, sBounds); // for strawSurface
      surface->assignSurfaceMaterial(surfaceMaterial);
      auto detElement = std::make_shared<MyDetectorElement>(std::make_shared<const Acts::Transform3>(stransform), surface, thickness);
      surface->assignDetectorElement(*detElement);
      detectorStore.push_back(std::move(detElement));
      vSurfaces.push_back(std::move(surface));
    }
    auto trafo = Acts::Transform3::Identity();
    trafo.translate(Acts::Vector3(0., 0., positions[i] * cm));
    const auto rBounds = std::make_shared<const Acts::RadialBounds>(rmin - 10, rmax + 10); // TODO check layer dimensions
    
    // creating custom surface array (surfaceArrayCreator.surfaceArrayOnDisc doesn't work)
    double tol = 1.;
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Closed> axisPhi(-M_PI, M_PI, 18); // TODO check binning
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Bound> axisR(rmin, rmax, 1);
    using SGL = Acts::SurfaceArray::SurfaceGridLookup<decltype(axisR), decltype(axisPhi)>;
    std::vector<Acts::AxisDirection> axisDirections = {Acts::AxisDirection::AxisR, Acts::AxisDirection::AxisPhi};
    auto repsurface = Acts::Surface::makeShared<Acts::DiscSurface>(trafo, rBounds);
    auto gridLookup = std::make_unique<SGL>(repsurface, tol, std::pair{axisR, axisPhi}, axisDirections);
    std::vector<const Acts::Surface*> surfacesRaw = unpackSmartPointers(vSurfaces);
    gridLookup->fill(gctx, surfacesRaw);
    auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(std::move(gridLookup),vSurfaces, trafo));

//    auto layer = Acts::DiscLayer::create(trafo, rBounds, std::move(surArray), thickness, nullptr, Acts::active);
    auto layer = Acts::DiscLayer::create(trafo, rBounds, std::move(surArray), 10.0, nullptr, Acts::active);

    for (auto& surface : vSurfaces) {
      auto mutableSurface = const_cast<Acts::Surface*>(surface.get());
      mutableSurface->associateLayer(*layer);
    }

    layerOrderVec.emplace_back(layer, layer->referencePosition(gctx, Acts::AxisDirection::AxisZ));
  }
  auto binning = std::make_unique<const Acts::BinUtility>(layerBinEdges, Acts::open, Acts::AxisDirection::AxisZ, Acts::Transform3::Identity());
  auto layArr = std::make_unique<Acts::BinnedArrayXD<Acts::LayerPtr>>(layerOrderVec, std::move(binning));
 
  auto boundsVol = std::make_shared<Acts::CuboidVolumeBounds>(2000._mm, 2000._mm, 3100._mm);
  auto trackVolume = std::make_shared<Acts::TrackingVolume>(Acts::Transform3::Identity(), boundsVol, volumeMaterial, std::move(layArr), nullptr, Acts::MutableTrackingVolumeVector{}, "Telescope");
  auto trackingGeometry = new Acts::TrackingGeometry(trackVolume);

  // Acts::ObjVisualization3D objVis;
  // const Acts::TrackingVolume& tgVolume = *(trackingGeometry->highestTrackingVolume());
  // Acts::GeometryView3D::drawTrackingVolume(objVis, tgVolume, gctx);

  // Check geometry
  bool print_surface_info = 1;
  if (print_surface_info) {
    const Acts::TrackingVolume *highestTrackingVolume = trackingGeometry->highestTrackingVolume();
    std::cout << "  volumeId = " << highestTrackingVolume->geometryId() << std::endl;
    const Acts::LayerArray *confinedLayers = highestTrackingVolume->confinedLayers();
    for (const auto &layer : confinedLayers->arrayObjects())  {
      bool isDisk = dynamic_cast<const Acts::DiscLayer*>(layer.get()) != nullptr;
      bool isPlane = dynamic_cast<const Acts::PlaneLayer*>(layer.get()) != nullptr;
      const char* shape = isDisk ? "disk" : (isPlane ? "plane" : "unknown");
      std::cout << "  layerId = " << layer->geometryId();
      printf("  thickness = %f type = %d shape = %s\n", layer->thickness(), layer->layerType(), shape);

      if (layer->layerType()==-1) {
        const Acts::NavigationLayer* nlayer = dynamic_cast<const Acts::NavigationLayer*>(layer.get());
        //printf(" navigation %p\n", nlayer);
      } else if (layer->layerType()==1 && isDisk) {
        const Acts::DiscLayer* player = dynamic_cast<const Acts::DiscLayer*>(layer.get());
        // auto pos = Acts::Vector3(920,4,0);
        // auto dir = Acts::Vector3(0,0,1);
        // auto v = layer->surfaceArray()->neighbors(pos, dir);
        // printf("    v.size()=%d\n",v.size());
        // for (auto& s : v) { std::cout << "    v.surfaceId = " << s->geometryId() << std::endl;  }
        // std::cout << std::endl;
      } else if (layer->layerType()==1 && isPlane) {
        const Acts::PlaneLayer* player = dynamic_cast<const Acts::PlaneLayer*>(layer.get());
        // printf("   plane %p\n", player);
      }
      if (!layer->surfaceArray())
        continue;
      // for (const auto &surface : layer->surfaceArray()->surfaces())
      // {
      //   std::cout << "    surfaceId = " << surface->geometryId();
      // }
      // std::cout << std::endl;
    } //for layers
  } // print_surface_info


  return trackingGeometry;
}

#endif