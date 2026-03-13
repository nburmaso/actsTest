#include <memory>

#include "Acts/Geometry/CylinderLayer.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include "TString.h"

#include "MyFtdDetector.h"

using namespace Acts::UnitConstants;
using namespace Acts::UnitLiterals;
using std::numbers::pi;

MyFtdDetector::MyFtdDetector()
{
  fFtdGeo = std::make_shared<MyFtdGeo>();
  fSecGeo = std::make_shared<BaseTpcSectorGeo>();

  auto vacuum = Acts::Material::fromMassDensity(1e+10_cm, 1e+10_cm, 1.0000, 1, 1e-10_g / 1_cm3);
  fVolumeMaterial = std::make_shared<Acts::HomogeneousVolumeMaterial>(vacuum);
}

std::shared_ptr<const Acts::TrackingGeometry> MyFtdDetector::GetTrackingGeometry(bool addFtd, bool addROC, bool addFlange, bool addTpc)
{
  if (fTrackingGeometry == nullptr) {
    CreateTrackingGeometry(addFtd, addROC, addFlange, addTpc);
  }
  return fTrackingGeometry;
}

std::shared_ptr<Acts::TrackingVolume> MyFtdDetector::CreateEndCap(bool isPos, bool addROC, bool addFlange)
{
  double sign = isPos ? 1 : -1; // positive/negative z

  std::unique_ptr<const Acts::LayerArray> layArr = nullptr;

  double zMinEndCap = fFtdGeo->GetZMinEndCap() * cm;
  double zMaxEndCap = fFtdGeo->GetZMaxEndCap() * cm;
  double rMax = fFtdGeo->GetRMax() * cm;

  if (addROC) {
    double rocThick = fFtdGeo->GetRocThick() * cm;
    double rMinRoc = fFtdGeo->GetRocRMin() * cm;
    double rMaxRoc = fFtdGeo->GetRocRMax() * cm;
    double zRoc = zMinEndCap + rocThick * 0.5 + 1.;
    int nSectors = fSecGeo->GetSectorCountHalf();

    Acts::LayerVector layVecEndCap;

    auto aluminium = Acts::Material::fromMassDensity(8.897_cm, 39.70_cm, 26.9815, 13, 2.699_g / 1_cm3);
    Acts::MaterialSlab matPropROC(aluminium, rocThick);
    const auto surfaceMatROC = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matPropROC);
    Acts::Translation3 trans(0, 0, zRoc);
    Acts::Transform3 trafo(trans);
    const auto rBoundsROC = std::make_shared<const Acts::RadialBounds>(rMinRoc, rMaxRoc);
    auto surfROC = Acts::Surface::makeShared<Acts::DiscSurface>(trafo, rMinRoc, rMaxRoc);
    surfROC->assignSurfaceMaterial(surfaceMatROC);
    auto surArray = std::make_unique<Acts::SurfaceArray>(surfROC);
    auto layer = Acts::DiscLayer::create(trafo, rBoundsROC, std::move(surArray), rocThick, nullptr, Acts::passive);
    surfROC->associateLayer(*layer);
    layVecEndCap.push_back(layer);

    if (addFlange) {
      double frameThick = fFtdGeo->GetFrameThick() * cm;
      double zFrame = zMinEndCap + rocThick + frameThick * 0.5 + 2.; // mm; todo: correct drift length????
      Acts::MaterialSlab matPropFrame(aluminium, frameThick);
      const auto surfaceMatFrame = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matPropFrame);
      std::vector<std::shared_ptr<const Acts::Surface>> vSurfaceFrame;

      Acts::Translation3 transFrameCircum1(0, 0, sign * zFrame);
      Acts::Translation3 transFrameCircum2(0, 0, sign * zFrame);
      Acts::Transform3 trafoFrameCircum1(transFrameCircum1);
      Acts::Transform3 trafoFrameCircum2(transFrameCircum2);

      const auto pBoundsFrameRadial = std::make_shared<const Acts::RectangleBounds>(fFtdGeo->GetFrameHalfXRadial() * cm,
                                                                                    fFtdGeo->GetFrameHalfYRadial() * cm);
      for (int iSector = 0; iSector < nSectors; iSector++) {
        Acts::Transform3 trafoFrameRadial;
        trafoFrameRadial.setIdentity();
        double phi = Acts::detail::radian_sym(fSecGeo->SectorAxisAngleRad(iSector));
        trafoFrameRadial.rotate(Eigen::AngleAxisd(pi / 12. + phi, Acts::Vector3(0, 0, 1)));
        trafoFrameRadial.translate(Acts::Vector3(fFtdGeo->GetFrameCXRadial0() * cm, 0, sign * zFrame));
        auto surfFrame = Acts::Surface::makeShared<Acts::PlaneSurface>(trafoFrameRadial, pBoundsFrameRadial);
        surfFrame->assignSurfaceMaterial(surfaceMatFrame);
        vSurfaceFrame.push_back(surfFrame);
        auto bounds1 = std::make_shared<const Acts::RadialBounds>(fFtdGeo->GetFrameRMin1() * cm,fFtdGeo->GetFrameRMax1() * cm, fSecGeo->GetSectorPhiRad() * 0.5 - 1e-10, phi);
        auto bounds2 = std::make_shared<const Acts::RadialBounds>(fFtdGeo->GetFrameRMin2() * cm,fFtdGeo->GetFrameRMax2() * cm, fSecGeo->GetSectorPhiRad() * 0.5 - 1e-10, phi);
        auto surfaceFrameCircum1 = Acts::Surface::makeShared<Acts::DiscSurface>(trafoFrameCircum1, bounds1);
        auto surfaceFrameCircum2 = Acts::Surface::makeShared<Acts::DiscSurface>(trafoFrameCircum2, bounds2);
        surfaceFrameCircum1->assignSurfaceMaterial(surfaceMatFrame);
        surfaceFrameCircum2->assignSurfaceMaterial(surfaceMatFrame);
        vSurfaceFrame.push_back(surfaceFrameCircum1);
        vSurfaceFrame.push_back(surfaceFrameCircum2);
      }

      Acts::SurfaceArrayCreator::Config sacConfig;
      Acts::SurfaceArrayCreator surfaceArrayCreator(sacConfig, Acts::getDefaultLogger("SAC-ECAP", Acts::Logging::INFO));
      const auto rBoundsFrame = std::make_shared<const Acts::RadialBounds>(fFtdGeo->GetFrameRMin1() * cm, fFtdGeo->GetFrameRMax2() * cm);
      Acts::Translation3 transFrame(0., 0., sign * zFrame);
      Acts::Transform3 trafoFrame(transFrame);
      auto layerFrame = Acts::DiscLayer::create(trafoFrame, rBoundsFrame, surfaceArrayCreator.surfaceArrayOnDisc(fGeoCtx, vSurfaceFrame, 3, 12), frameThick, nullptr, Acts::passive);

      for (auto& surface : vSurfaceFrame) {
        auto mutableSurface = const_cast<Acts::Surface*>(surface.get());
        mutableSurface->associateLayer(*layerFrame);
      }
      layVecEndCap.push_back(layerFrame);
    }

    // Create layer array
    Acts::LayerArrayCreator::Config lacConfig;
    Acts::LayerArrayCreator layArrCreator(lacConfig, Acts::getDefaultLogger("LAC-ECAP", Acts::Logging::INFO));
    double zMin, zMax;
    if (sign > 0) {
      zMin = zMinEndCap;
      zMax = zMaxEndCap;
    } else {
      zMin = -zMaxEndCap;
      zMax = -zMinEndCap;
    }
    layArr = std::unique_ptr<const Acts::LayerArray>(layArrCreator.layerArray(fGeoCtx, layVecEndCap, zMin, zMax, Acts::BinningType::arbitrary, Acts::AxisDirection::AxisZ));
  }

  Acts::Transform3 trafoEndCap = Acts::Transform3::Identity();
  trafoEndCap.translate(Acts::Vector3(0., 0., sign * (zMaxEndCap + zMinEndCap) / 2.));
  auto endCapBounds = std::make_shared<Acts::CylinderVolumeBounds>(0, rMax, (zMaxEndCap - zMinEndCap) / 2.);
  auto volName = isPos ? "PosEndCap" : "NegEndCap";
  auto endCapVolume = std::make_shared<Acts::TrackingVolume>(trafoEndCap, endCapBounds, fVolumeMaterial,
                                                             std::move(layArr), nullptr,
                                                             Acts::MutableTrackingVolumeVector{}, volName);

  return endCapVolume;
}

std::shared_ptr<Acts::TrackingVolume> MyFtdDetector::CreateFTD(bool isPos, bool addLayers)
{
  double sign = isPos ? 1 : -1; // positive/negative z
  auto volName = isPos ? "PosFTD" : "NegFTD";

  double rMax = fFtdGeo->GetRMax() * cm;
  double zMinFtd = fFtdGeo->GetFtdZMin() * cm;
  double zMaxFtd = fFtdGeo->GetFtdZMax() * cm;
  Acts::Transform3 trafoFtd = Acts::Transform3::Identity();
  trafoFtd.translate(Acts::Vector3(0., 0., sign * (zMaxFtd + zMinFtd) / 2.));
  auto ftdBounds = std::make_shared<Acts::CylinderVolumeBounds>(0, 20000, (zMaxFtd - zMinFtd) / 2.); // rmin, rmax, halfz

  std::unique_ptr<const Acts::LayerArray> layArr = nullptr;

  if (addLayers) { // TODO adjust for negative sign
    std::vector<float> layerBinEdges;
    std::vector<std::pair<std::shared_ptr<const Acts::Layer>, Acts::Vector3>> layerOrderVec;
    const auto& positions = fFtdGeo->GetLayerPositions();
    auto silicon = Acts::Material::fromMassDensity(9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);

    // add first navigation layer
    Acts::Transform3 tr1 = Acts::Transform3::Identity();
    double z1 = zMinFtd;
    double z2 = (positions[0]*cm - fFtdGeo->GetLayerThickness(0)*cm/2.);
    tr1.translate(Acts::Vector3(0, 0, (z2 + z1)/2.));
    const auto pBounds = std::make_shared<const Acts::RectangleBounds>(rMax, rMax);
    auto layer = Acts::PlaneLayer::create(tr1, pBounds, nullptr, z2 - z1, nullptr, Acts::navigation);
    layerOrderVec.emplace_back(layer, layer->referencePosition(fGeoCtx, Acts::AxisDirection::AxisZ));
    layerBinEdges.push_back(z1);
    layerBinEdges.push_back(z2);

    for (int i = 0; i < positions.size(); ++i) {
      printf("Adding layer %d type %d\n",i, fFtdGeo->GetLayerType(i));
      double rmin = fFtdGeo->GetLayerRMin(i) * cm;
      double rmax = fFtdGeo->GetLayerRMax(i) * cm;
      layerBinEdges.push_back(i < positions.size()-1 ? (positions[i] + positions[i+1])/2.*cm : zMaxFtd); // TODO sign
      double thickness = fFtdGeo->GetLayerThickness(i) * cm;
      // create surface material
      auto matProp = Acts::MaterialSlab(silicon, thickness);
      auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);

      if (fFtdGeo->GetLayerType(i)==2) {
        const auto pBounds = std::make_shared<const Acts::RectangleBounds>(rmax, rmax);
        Acts::Transform3 trafo = Acts::Transform3::Identity();
        trafo.translate(Acts::Vector3(0, 0, sign*positions[i]*cm));
        auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafo, pBounds);
        surface->assignSurfaceMaterial(std::move(surfaceMaterial));
        auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(surface));
        auto layer = Acts::PlaneLayer::create(trafo, pBounds, std::move(surArray), thickness, nullptr, Acts::active);
        surface->associateLayer(*layer.get());
        layerOrderVec.emplace_back(layer, layer->referencePosition(fGeoCtx, Acts::AxisDirection::AxisZ));
        auto detElement = std::make_shared<MyFtdDetectorElement>(std::make_shared<const Acts::Transform3>(trafo), surface, thickness);
        detElement->setName(Form("%s%d", volName, i));
        detElement->setLayer(i);
        surface->assignDetectorElement(*detElement);
        fDetectorStoreFtd.push_back(std::move(detElement));
        continue;
      }
      // TODO: use ftdGeo getters for angles, rc etc.
      double angleRot = 0;
      if (fFtdGeo->GetLayerType(i)==4) angleRot = 0;
      if (fFtdGeo->GetLayerType(i)==5) angleRot = +fFtdGeo->GetTubeIncl();
      if (fFtdGeo->GetLayerType(i)==6) angleRot = -fFtdGeo->GetTubeIncl();
      double rc = 0.5 * (rmax + rmin);
      double hl = 0.5 * (rmax - rmin) / cos(angleRot);
      double hw = fFtdGeo->GetLayerStripWidth(i) * cm;
      const auto sBounds = std::make_shared<const Acts::RectangleBounds>(hw, hl);
      int nTubes = fFtdGeo->GetLayerNumberOfTubes(i);
      std::vector<std::shared_ptr<const Acts::Surface>> vSurfaces;
      for (int iTube = 0; iTube < nTubes; iTube++) {
        double phi = fFtdGeo->GetLayerAngle(i) + 2*pi/nTubes*(iTube + 0.5);
        double cosp = cos(phi);
        double sinp = sin(phi);
        auto stransform = Acts::Transform3::Identity();
        stransform.translate(Acts::Vector3(rc*cosp, rc*sinp, sign*positions[i] * cm));
        stransform.rotate(Eigen::AngleAxisd(pi/2 + phi + angleRot, Acts::Vector3(0, 0, 1)));
        auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(stransform, sBounds);
        surface->assignSurfaceMaterial(surfaceMaterial);
        auto detElement = std::make_shared<MyFtdDetectorElement>(std::make_shared<const Acts::Transform3>(stransform), surface, thickness);
        detElement->setName(Form("%s%d_%d", volName, i, iTube));
        detElement->setLayer(i);
        surface->assignDetectorElement(*detElement);
        fDetectorStoreFtd.push_back(std::move(detElement));
        vSurfaces.push_back(std::move(surface));
      }
      auto trafo = Acts::Transform3::Identity();
      trafo.translate(Acts::Vector3(0., 0., sign*positions[i] * cm));
      const auto rBounds = std::make_shared<const Acts::RadialBounds>(rmin - 10, rmax + 10); // TODO check layer dimensions

      // creating custom surface array (surfaceArrayCreator.surfaceArrayOnDisc doesn't work)
      double tol = 1.;
      Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Closed> axisPhi(-pi, pi, 18); // TODO check binning
      Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Bound> axisR(rmin, rmax, 1);
      using SGL = Acts::SurfaceArray::SurfaceGridLookup<decltype(axisR), decltype(axisPhi)>;
      std::vector<Acts::AxisDirection> axisDirections = {Acts::AxisDirection::AxisR, Acts::AxisDirection::AxisPhi};
      auto repsurface = Acts::Surface::makeShared<Acts::DiscSurface>(trafo, rBounds);
      auto gridLookup = std::make_unique<SGL>(repsurface, tol, std::pair{axisR, axisPhi}, axisDirections);
      std::vector<const Acts::Surface*> surfacesRaw = unpackSmartPointers(vSurfaces);
      gridLookup->fill(fGeoCtx, surfacesRaw);
      auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(std::move(gridLookup),vSurfaces, trafo));
      auto layer = Acts::DiscLayer::create(trafo, rBounds, std::move(surArray), thickness, nullptr, Acts::active);

      for (auto& surface : vSurfaces) {
        auto mutableSurface = const_cast<Acts::Surface*>(surface.get());
        mutableSurface->associateLayer(*layer);
      }
      layerOrderVec.emplace_back(layer, layer->referencePosition(fGeoCtx, Acts::AxisDirection::AxisZ));
    }
    // sort layerBinEdges for negative z?
    auto binning = std::make_unique<const Acts::BinUtility>(layerBinEdges, Acts::open, Acts::AxisDirection::AxisZ, Acts::Transform3::Identity());
    layArr = std::make_unique<Acts::BinnedArrayXD<Acts::LayerPtr>>(layerOrderVec, std::move(binning));
  }

  auto ftdVolume = std::make_shared<Acts::TrackingVolume>(trafoFtd, std::move(ftdBounds), fVolumeMaterial, std::move(layArr), nullptr, Acts::MutableTrackingVolumeVector{}, volName);
  return ftdVolume;
}

void MyFtdDetector::CreateTrackingGeometry(bool addFtd, bool addROC, bool addFlange, bool addTpc)
{
  const auto& rowCount = fSecGeo->GetRowCount();
  const auto& padHeight = fSecGeo->GetPadHeight();
  const auto& padWidth = fSecGeo->GetPadWidth();
  const std::vector<int>& halfPadCount = fSecGeo->GetPadCount();
  int nPadRows = halfPadCount.size();

  double mpdHalfZ = fFtdGeo->GetZMaxEndCap() * cm;
  double mpdVolR = fFtdGeo->GetRMax() * cm;

  double barrelHalfZ = fFtdGeo->GetZMinEndCap() * cm;
  double innerTpcVolumeR = fSecGeo->GetYPadAreaLowerEdge() * cm - padHeight[0] * cm;
  double outerTpcVolumeR = mpdVolR;

  double pipeVolumeR = innerTpcVolumeR;
  double pipeR = fFtdGeo->GetPipeR() * cm;
  double pipeThick = fFtdGeo->GetPipeThick() * cm;

  int nSectors = fSecGeo->GetSectorCountHalf();
  double tpcSecHalfPhi = 0.5 * fSecGeo->GetSectorPhiRad();

  // Create materials
  // todo: fill pipe with air
  // auto air = Acts::Material::fromMassDensity(3.039e+4_cm, 7.477e+4_cm, 28.97, 28.97 * 0.49919, 1.205e-3_g / 1_cm3); // dry, 1 atm, 20C, todo: fetch from geant?
  auto argon = Acts::Material::fromMassDensity(1.176e+4_cm, 7.204e+4_cm, 39.948, 18, 1.662e-3_g / 1_cm3);

  auto argonMatPropInner = Acts::MaterialSlab(argon, padHeight[0] * cm);
  auto argonMatInner = std::make_shared<Acts::HomogeneousSurfaceMaterial>(argonMatPropInner);
  auto argonMatPropOuter = Acts::MaterialSlab(argon, padHeight[1] * cm);
  auto argonMatOuter = std::make_shared<Acts::HomogeneousSurfaceMaterial>(argonMatPropOuter);
  auto argonVolumeMaterial = std::make_shared<Acts::HomogeneousVolumeMaterial>(argon);

  // create endcap and FTD volumes
  // using function to avoid duplication for two z directions
  std::shared_ptr<Acts::TrackingVolume> posEndCapVolume = CreateEndCap(true, addROC, addFlange);
  std::shared_ptr<Acts::TrackingVolume> posFtdVolume = CreateFTD(true, addFtd);

  // create pipe volume
  // todo: fetch geometry and material parameters from pipe root file
  auto beryllium = Acts::Material::fromMassDensity(35.28_cm, 42.10_cm, 9.0122, 4, 1.848_g / 1_cm3);
  auto surfaceMaterialPipe = std::make_shared<Acts::HomogeneousSurfaceMaterial>(Acts::MaterialSlab(beryllium, pipeThick));

  const auto pipeSurfBounds = std::make_shared<const Acts::CylinderBounds>(pipeR, barrelHalfZ);
  auto pipeSurface = Acts::Surface::makeShared<Acts::CylinderSurface>(Acts::Transform3::Identity(), pipeSurfBounds);
  pipeSurface->assignSurfaceMaterial(surfaceMaterialPipe);
  auto pipeSurfArr = std::make_unique<Acts::SurfaceArray>(pipeSurface);
  const auto pipeLayerBounds = std::make_shared<const Acts::CylinderBounds>(pipeR, barrelHalfZ);
  auto layerPipe = Acts::CylinderLayer::create(Acts::Transform3::Identity(), pipeLayerBounds, std::move(pipeSurfArr), pipeThick, nullptr, Acts::passive);
  Acts::LayerVector pipeLayVec{layerPipe};
  pipeSurface->associateLayer(*layerPipe);
  Acts::LayerArrayCreator pipeLayArrCreator({}, Acts::getDefaultLogger("LAC-PIPE", Acts::Logging::INFO));
  std::unique_ptr<const Acts::LayerArray> pipeLayArr(pipeLayArrCreator.layerArray(fGeoCtx, pipeLayVec, 0., pipeVolumeR, Acts::BinningType::arbitrary, Acts::AxisDirection::AxisR));
  auto pipeVolumeBounds = std::make_shared<Acts::CylinderVolumeBounds>(0., pipeVolumeR, barrelHalfZ);
  auto pipeVolume = std::make_shared<Acts::TrackingVolume>(Acts::Transform3::Identity(), pipeVolumeBounds, fVolumeMaterial, std::move(pipeLayArr), nullptr, Acts::MutableTrackingVolumeVector{}, "PIPE");

  // below: create TPC volume

  // containers that are used if TPC sectors need to be created
  Acts::TrackingVolumeVector tpcSecVolVec;
  std::shared_ptr<const Acts::TrackingVolumeArray> tpcSecVolArr = nullptr;
  std::shared_ptr<const Acts::BinnedArrayXD<Acts::TrackingVolumePtr>> binnedSectorArray = nullptr;
  const std::vector<int> sectorOrder{6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5}; // filling TPC with sectors from -pi to +pi -> following ACTS conventions

  // create sector volumes if needed -- fill TPC with material otherwise
  if (addTpc) {
    for (int sec : sectorOrder) {
      double tpcSecAvgPhi = Acts::detail::radian_sym(fSecGeo->SectorAxisAngleRad(sec));
      Acts::Transform3 trafoTpcSec = Acts::Transform3::Identity();
      trafoTpcSec.rotate(Eigen::AngleAxisd(tpcSecAvgPhi, Acts::Vector3(0, 0, 1)));
      double yShift = fSecGeo->GetYPadAreaLowerEdge();
      std::vector<float> layerBinEdges;
      std::vector<std::pair<std::shared_ptr<const Acts::Layer>, Acts::Vector3>> tpcLayOrderVec;
      for (int i = 0; i < nPadRows; ++i) {
        double rowCenterY = (yShift + fSecGeo->PadRowCenter2Local(0, i).Y()) * cm;
        bool isInner = i < rowCount[0];
        double padW = padWidth[isInner ? 0 : 1] * cm;
        double padH = padHeight[isInner ? 0 : 1] * cm;
        double rowHalfWidth = halfPadCount[i] * padW;
        const auto layerBounds = std::make_shared<const Acts::RectangleBounds>(rowHalfWidth, barrelHalfZ);
        auto trafoPadRow = Acts::Transform3::Identity();
        trafoPadRow.rotate(Eigen::AngleAxisd(pi / 2., Acts::Vector3(1, 0, 0)));
        trafoPadRow.rotate(Eigen::AngleAxisd(pi / 2., Acts::Vector3(0, 1, 0)));
        trafoPadRow.rotate(Eigen::AngleAxisd(tpcSecAvgPhi, Acts::Vector3(0, 1, 0)));
        trafoPadRow.translate(Acts::Vector3(0., 0., rowCenterY));

        auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafoPadRow, layerBounds);
        surface->assignSurfaceMaterial(isInner ? argonMatInner : argonMatOuter);

        auto detElement = std::make_shared<MyTpcDetectorElement>(std::make_shared<const Acts::Transform3>(trafoPadRow), surface, padH);
        detElement->setName(Form("TPC_sec%d_row%d", sec, i));
        detElement->setSector(sec);
        detElement->setPadRow(i);
        surface->assignDetectorElement(*detElement);
        fDetectorStoreTpc.push_back(std::move(detElement));

        auto surfArr = std::make_unique<Acts::SurfaceArray>(surface);
        auto layer = Acts::PlaneLayer::create(trafoPadRow, layerBounds, std::move(surfArr), padH, nullptr, Acts::active);
        surface->associateLayer(*layer);

        tpcLayOrderVec.emplace_back(layer, layer->referencePosition(fGeoCtx, Acts::AxisDirection::AxisR));
        auto r = layer->referencePositionValue(fGeoCtx, Acts::AxisDirection::AxisR);
        if (i == 0) layerBinEdges.push_back(r - padH / 2.);
        layerBinEdges.push_back(r + padH / 2.);
      }
      auto binning = std::make_unique<const Acts::BinUtility>(layerBinEdges, Acts::open, Acts::AxisDirection::AxisR, trafoTpcSec);
      auto layArr = std::make_unique<Acts::BinnedArrayXD<Acts::LayerPtr>>(tpcLayOrderVec, std::move(binning));
      auto tpcSectorBounds = std::make_shared<Acts::CylinderVolumeBounds>(innerTpcVolumeR, outerTpcVolumeR, barrelHalfZ, tpcSecHalfPhi, 0.);
      auto tpcSectorVolume = std::make_shared<Acts::TrackingVolume>(trafoTpcSec, tpcSectorBounds, fVolumeMaterial,
                                                                    std::move(layArr), nullptr, Acts::MutableTrackingVolumeVector{},
                                                                    Form("TPC%d", sec));
      tpcSecVolVec.push_back(tpcSectorVolume);
    }

    // glue TPC sector volumes with each other
    for (int i = 0; i < nSectors; ++i) {
      auto tpcSecVol1 = std::const_pointer_cast<Acts::TrackingVolume>(tpcSecVolVec[i]);
      auto tpcSecVol2 = std::const_pointer_cast<Acts::TrackingVolume>(tpcSecVolVec[(i + 1) % nSectors]);
      tpcSecVol2->glueTrackingVolume(fGeoCtx, Acts::BoundarySurfaceFace::tubeSectorNegativePhi, tpcSecVol1.get(), Acts::BoundarySurfaceFace::tubeSectorPositivePhi);
    }

    // create sector volume array
    Acts::TrackingVolumeArrayCreator::Config tpcVolArrayCreatorCfg;
    Acts::TrackingVolumeArrayCreator tpcVolArrayCreator(tpcVolArrayCreatorCfg);
    tpcSecVolArr = tpcVolArrayCreator.trackingVolumeArray(fGeoCtx, tpcSecVolVec, Acts::AxisDirection::AxisPhi);
  }

  // create TPC volume -- contains TPC sector volumes if created
  auto tpcVolBounds = std::make_shared<Acts::CylinderVolumeBounds>(innerTpcVolumeR, outerTpcVolumeR, barrelHalfZ);
  auto tpcMaterial = tpcSecVolArr == nullptr ? argonVolumeMaterial : fVolumeMaterial; // fill with material if no sectors were created
  auto tpcVolume = std::make_shared<Acts::TrackingVolume>(Acts::Transform3::Identity(), tpcVolBounds, tpcMaterial, nullptr,
                                                          std::move(tpcSecVolArr), Acts::MutableTrackingVolumeVector{}, "TPC");

  // join pipe and tpc
  tpcVolume->glueTrackingVolume(fGeoCtx, Acts::BoundarySurfaceFace::tubeInnerCover, pipeVolume.get(), Acts::BoundarySurfaceFace::cylinderCover);
  tpcVolume->glueTrackingVolume(fGeoCtx, Acts::BoundarySurfaceFace::positiveFaceXY, posEndCapVolume.get(), Acts::BoundarySurfaceFace::negativeFaceXY);
  pipeVolume->glueTrackingVolume(fGeoCtx, Acts::BoundarySurfaceFace::positiveFaceXY, posEndCapVolume.get(), Acts::BoundarySurfaceFace::negativeFaceXY);

  // attach TPC volumes to the pipe and endcap outer surface
  if (addTpc) {
    // create phi binning for tpc sector volumes
    std::vector<std::pair<Acts::TrackingVolumePtr, Acts::Vector3>> vTpcSecVolumesBinning;
    std::vector<float> tpcSecBounds;
    for (int i = 0; i < nSectors; ++i) {
      const auto& tpcSecVol = tpcSecVolVec[i];
      int sec = sectorOrder[i];
      double tpcSecAvgPhi = Acts::detail::radian_sym(fSecGeo->SectorAxisAngleRad(sec));
      vTpcSecVolumesBinning.emplace_back(tpcSecVol, Acts::Vector3(std::cos(tpcSecAvgPhi), std::sin(tpcSecAvgPhi), 0.));
      if (i == 0) tpcSecBounds.push_back(tpcSecAvgPhi - tpcSecHalfPhi);
      tpcSecBounds.push_back(tpcSecAvgPhi + tpcSecHalfPhi);
    }

    auto binning = std::make_unique<const Acts::BinUtility>(tpcSecBounds, Acts::BinningOption::closed, Acts::AxisDirection::AxisPhi);
    binnedSectorArray = std::make_shared<const Acts::BinnedArrayXD<Acts::TrackingVolumePtr>>(vTpcSecVolumesBinning, std::move(binning));

    auto pipeBoundSurface = pipeVolume->boundarySurfaces().at(Acts::BoundarySurfaceFace::cylinderCover);
    auto mutablePipeBoundSurface = std::const_pointer_cast<Acts::BoundarySurfaceT<Acts::TrackingVolume>>(pipeBoundSurface);
    mutablePipeBoundSurface->attachVolumeArray(binnedSectorArray, Acts::Direction::Positive());

    auto posEndCapSurface = posEndCapVolume->boundarySurfaces().at(Acts::BoundarySurfaceFace::negativeFaceXY);
    auto mutablePosEndCapSurface = std::const_pointer_cast<Acts::BoundarySurfaceT<Acts::TrackingVolume>>(posEndCapSurface);
    mutablePosEndCapSurface->attachVolumeArray(binnedSectorArray, Acts::Direction::Negative());

    for (const auto& tpcSecVol : tpcSecVolVec) {
      auto mutableSectorVol = std::const_pointer_cast<Acts::TrackingVolume>(tpcSecVol);
      auto sectorSurfaceInnerCover = mutableSectorVol->boundarySurfaces().at(Acts::BoundarySurfaceFace::tubeInnerCover);
      auto mutableSectorSurfaceInnerCover = std::const_pointer_cast<Acts::BoundarySurfaceT<Acts::TrackingVolume>>(sectorSurfaceInnerCover);
      auto sectorSurfaceFace = mutableSectorVol->boundarySurfaces().at(Acts::BoundarySurfaceFace::positiveFaceXY);
      auto mutableSectorSurfaceFace = std::const_pointer_cast<Acts::BoundarySurfaceT<Acts::TrackingVolume>>(sectorSurfaceFace);
      mutableSectorSurfaceInnerCover->attachVolume(pipeVolume.get(), Acts::Direction::Negative()); // sector->pipe
      mutableSectorSurfaceFace->attachVolume(posEndCapVolume.get(), Acts::Direction::Positive()); // sector->endcap
    }
  }

  // put pipe and TPC into barrel volume
  Acts::TrackingVolumeArrayCreator::Config barrelArrayCreatorCfg;
  Acts::TrackingVolumeArrayCreator barrelArrayCreator(barrelArrayCreatorCfg);
  Acts::TrackingVolumeVector vBarrelVolumes = {pipeVolume, tpcVolume};
  auto barrelArray = barrelArrayCreator.trackingVolumeArray(fGeoCtx, vBarrelVolumes, Acts::AxisDirection::AxisR);
  auto barrelBounds = std::make_shared<Acts::CylinderVolumeBounds>(0., mpdVolR, barrelHalfZ);
  auto barrelVolume = std::make_shared<Acts::TrackingVolume>(Acts::Transform3::Identity(), barrelBounds, fVolumeMaterial,
                                                             nullptr, std::move(barrelArray),
                                                             Acts::MutableTrackingVolumeVector{}, "BARREL");

  // attach TPC sector volumes to Barrel surface
  if (addTpc) {
    auto barrelSurface = barrelVolume->boundarySurfaces().at(Acts::BoundarySurfaceFace::positiveFaceXY);
    auto mutableBarrelSurface = std::const_pointer_cast<Acts::BoundarySurfaceT<Acts::TrackingVolume>>(barrelSurface);
    mutableBarrelSurface->attachVolumeArray(binnedSectorArray, Acts::Direction::Negative());
  }

  // glue: Barrel->Endcap->FTD
  barrelVolume->glueTrackingVolume(fGeoCtx, Acts::BoundarySurfaceFace::positiveFaceXY, posEndCapVolume.get(), Acts::BoundarySurfaceFace::negativeFaceXY);
  posEndCapVolume->glueTrackingVolume(fGeoCtx, Acts::BoundarySurfaceFace::positiveFaceXY, posFtdVolume.get(), Acts::BoundarySurfaceFace::negativeFaceXY);

  // create mother tracking volume
  Acts::TrackingVolumeVector vMpdVolumes = {barrelVolume, posEndCapVolume, posFtdVolume};
  Acts::TrackingVolumeArrayCreator::Config mpdArrayCreatorCfg;
  Acts::TrackingVolumeArrayCreator mpdArrayCreator(mpdArrayCreatorCfg);
  auto mpdVolumeArray = mpdArrayCreator.trackingVolumeArray(fGeoCtx, vMpdVolumes, Acts::AxisDirection::AxisZ);
  auto mpdVolumeBounds = std::make_shared<Acts::CylinderVolumeBounds>(0, mpdVolR, mpdHalfZ);
  auto mpdVolume = std::make_shared<Acts::TrackingVolume>(Acts::Transform3::Identity(), mpdVolumeBounds, fVolumeMaterial,
                                                          nullptr, std::move(mpdVolumeArray),
                                                          Acts::MutableTrackingVolumeVector{}, "MPD");

  fTrackingGeometry = std::make_shared<const Acts::TrackingGeometry>(mpdVolume);

  // Acts::ObjVisualization3D objVis;
  // const Acts::TrackingVolume& tgVolume = *(fTrackingGeometry->highestTrackingVolume());
  // Acts::GeometryView3D::drawTrackingVolume(objVis, tgVolume, fGeoCtx);

  std::string momVolName = "";

  auto surfVisitor = [&momVolName, this](const Acts::Surface* srf) {
    std::cout  << srf->geometryId() << " ";
  };

  // go through volumes, print info, collect volume names and ACTS geometry IDs
  std::function<void(std::vector<std::shared_ptr<const Acts::TrackingVolume>>&, int)> volVisitor;
  volVisitor = [&momVolName, &volVisitor, &surfVisitor, this] (std::vector<std::shared_ptr<const Acts::TrackingVolume>>& volArray, int recLevel) {
    recLevel++;
    if (recLevel > 5)
      return;
    for (const auto& v : volArray) {
      std::string vName = v->volumeName();
      momVolName = vName;
      fVolNames[v->geometryId().volume()] = vName;
      std::cout << std::string(recLevel, ' ') << "volume" << recLevel << ": " << vName.data() << "  " << v->geometryId();
      for (const auto &boundary : v->boundarySurfaces()) {
        std::cout << "  " << boundary->surfaceRepresentation().geometryId();
      }
      std::cout << std::endl;
      if (momVolName.find("TPC0") != std::string::npos) {
        std::cout << std::string(recLevel, ' ');
        v->visitSurfaces(surfVisitor);
        std::cout << std::endl;
      }
      auto conf = v->confinedVolumes();
      if (conf != nullptr) {
        auto inVolumes = conf->arrayObjects();
        volVisitor(inVolumes, recLevel);
      }
    }
  };

  int recLevel = 0;
  auto volume1 = fTrackingGeometry->highestTrackingVolume()->confinedVolumes()->arrayObjects();
  volVisitor(volume1, recLevel);

  CollectTpcSecRowGeoId(fTrackingGeometry);
  CollectFtdLayerGeoId(fTrackingGeometry);
}
