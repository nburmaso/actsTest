#pragma once

#include <map>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"

#include "MyFtdGeo.h"
#include "BaseTpcSectorGeo.h"

class MyDetectorElement : public Acts::DetectorElementBase {
public:
  MyDetectorElement(std::shared_ptr<const Acts::Transform3> transform, std::shared_ptr<Acts::Surface> surface,
                    double thickness)
    : Acts::DetectorElementBase(), fTransform(transform), fSurface(surface), fThickness(thickness) {
  }

  const Acts::Transform3 &transform(const Acts::GeometryContext &gctx) const { return *fTransform; }

  const Acts::Surface &surface() const { return *fSurface; }

  Acts::Surface &surface() { return *fSurface; }

  double thickness() const { return fThickness; }

  const std::string &name() const { return fName; }

  void setName(const std::string &name) { fName = name; }

protected:
  std::shared_ptr<const Acts::Transform3> fTransform = nullptr;
  std::shared_ptr<Acts::Surface> fSurface = nullptr;
  double fThickness = 0.;
  std::string fName{""};
};

class MyFtdDetectorElement : public MyDetectorElement {
public:
  MyFtdDetectorElement(std::shared_ptr<const Acts::Transform3> transform, std::shared_ptr<Acts::Surface> surface,
                       double thickness)
    : MyDetectorElement(transform, surface, thickness) {
  }

  int layer() const { return fLayer; }

  void setLayer(int layer) { fLayer = layer; }

private:
  int fLayer{};
};

class MyTpcDetectorElement : public MyDetectorElement {
public:
  MyTpcDetectorElement(std::shared_ptr<const Acts::Transform3> transform, std::shared_ptr<Acts::Surface> surface,
                       double thickness)
    : MyDetectorElement(transform, surface, thickness) {
  }

  int sector() const { return fSector; }

  void setSector(int sector) { fSector = sector; }

  int padRow() const { return fPadRow; }

  void setPadRow(int padrow) { fPadRow = padrow; }

private:
  int fSector{};
  int fPadRow{};
};

class MyFtdDetector {

public:
  MyFtdDetector();

  ~MyFtdDetector() = default;

  std::shared_ptr<const Acts::TrackingGeometry>
  GetTrackingGeometry(bool addFtd = true, bool addROC = true, bool addFlange = true, bool addTpc = true);

  std::shared_ptr<MyFtdGeo> FtdGeo() { return fFtdGeo; }

  std::shared_ptr<BaseTpcSectorGeo> TpcGeo() { return fSecGeo; }

  const Acts::GeometryContext &GeoContext() { return fGeoCtx; }

  std::string VolumeName(int volumeId) {
    std::string res = "none";
    if (volumeId == fCacheVolId) {
      res = fCacheVolName;
    }
    auto it = fVolNames.find(volumeId);
    if (it != fVolNames.end()) {
      res = it->second;
      fCacheVolName = res;
    }
    return res;
  }

  Acts::GeometryIdentifier FtdLayerToGeoId(int layer) {
    if (layer == fCacheLay) {
      return fCacheFtdGeoId;
    }
    auto it = fFtdLayerToGeoId.find(layer);
    if (it != fFtdLayerToGeoId.end()) {
      fCacheFtdGeoId = it->second;
      return fCacheFtdGeoId;
    }
    Acts::GeometryIdentifier geoId;
    return geoId.withSensitive(-1);
  }

  int GeoIdToFtdLayer(Acts::GeometryIdentifier id) {
    static Acts::GeometryIdentifier chId{};
    static int chLayer = -1;
    if (id == chId) {
      return chLayer;
    }
    auto it = fGeoIdToFtdLayer.find(id);
    if (it != fGeoIdToFtdLayer.end()) {
      chLayer = it->second;
      chId = id;
    }
    return chLayer;
  }

private:
  // geometry context
  Acts::GeometryContext fGeoCtx;

  // materials
  std::shared_ptr<Acts::HomogeneousVolumeMaterial> fVolumeMaterial{nullptr};

  std::vector<std::shared_ptr<MyFtdDetectorElement>> fDetectorStoreFtd;
  std::vector<std::shared_ptr<MyTpcDetectorElement>> fDetectorStoreTpc;

  std::string fCacheVolName{"dummy"};
  int fCacheSec{-1};
  int fCacheRow{-1};
  int fCacheLay{-1};
  Acts::GeometryIdentifier fCacheTpcGeoId{};
  Acts::GeometryIdentifier fCacheFtdGeoId{};
  int fCacheVolId{-1};
  std::unordered_map<int, std::string> fVolNames{};
  std::map<std::pair<int, int>, Acts::GeometryIdentifier> fTpcSecRowToGeoId{};
  std::map<int, Acts::GeometryIdentifier> fFtdLayerToGeoId{};
  std::unordered_map<Acts::GeometryIdentifier, int> fGeoIdToFtdLayer{};
  std::unordered_map<Acts::GeometryIdentifier, std::pair<int, int>> fGeoIdToTpcSecRow{};

  std::shared_ptr<MyFtdGeo> fFtdGeo{nullptr};
  std::shared_ptr<BaseTpcSectorGeo> fSecGeo{nullptr};

  std::shared_ptr<const Acts::TrackingGeometry> fTrackingGeometry{nullptr};

  std::shared_ptr<Acts::TrackingVolume> CreateEndCap(bool isPos, bool addROC = true, bool addFlange = true);

  std::shared_ptr<Acts::TrackingVolume> CreateFTD(bool isPos, bool addLayers = true);

  void CreateTrackingGeometry(bool addFtd = true, bool addROC = true, bool addFlange = true, bool addTpc = true);

  void CollectTpcSecRowGeoId(const std::shared_ptr<const Acts::TrackingGeometry>& geo) {
    fTpcSecRowToGeoId.clear();
    fGeoIdToTpcSecRow.clear();
    const auto& surfaceByIdentifier = geo->geoIdSurfaceMap();
    for (auto surfId : surfaceByIdentifier) {
      const auto* detEl = dynamic_cast<const MyDetectorElement*>(surfId.second->associatedDetectorElement());
      if (detEl == nullptr)
        continue;
      if (detEl->name().find("TPC") == std::string::npos)
        continue;
      const auto* detElTpc = dynamic_cast<const MyTpcDetectorElement*>(detEl);
      int sector = detElTpc->sector();
      int row = detElTpc->padRow();
      auto id = surfId.first;
      fTpcSecRowToGeoId[{sector, row}] = id;
      fGeoIdToTpcSecRow[id] = std::make_pair(sector, row);
    }
  }

  void CollectFtdLayerGeoId(const std::shared_ptr<const Acts::TrackingGeometry>& geo) {
    fFtdLayerToGeoId.clear();
    fGeoIdToFtdLayer.clear();
    const auto& surfaceByIdentifier = geo->geoIdSurfaceMap();
    for (auto surfId : surfaceByIdentifier) {
      const auto* detEl = dynamic_cast<const MyDetectorElement*>(surfId.second->associatedDetectorElement());
      if (detEl == nullptr)
        continue;
      if (detEl->name().find("FTD") == std::string::npos)
        continue;
      const auto* detElFtd = dynamic_cast<const MyFtdDetectorElement*>(detEl);
      int layer = detElFtd->layer();
      auto id = surfId.first;
      fFtdLayerToGeoId[layer] = id;
      fGeoIdToFtdLayer[id] = layer;
    }
  }
};