// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "MySpacePointMaker.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderOptions.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
//#include "ActsExamples/Utilities/GroupBy.hpp"

#include <algorithm>
#include <functional>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <list>
#include <map>

#include "TMath.h"
//#include "tracker_config.h"
#include "Math/Functor.h"
#include "Math/RootFinder.h"

// helper namespace
namespace {
  struct FunctorFk {
    int n;
    double tx;
    double ty;
    std::vector<double>* z{nullptr};
    std::vector<double>* s{nullptr};
    std::vector<double>* c{nullptr};
    std::vector<double>* g{nullptr};
    std::vector<double>* sk{nullptr};
    std::vector<double>* ck{nullptr};
    double operator() (double kkk) {
      double ss = 0;
      double sc = 0;
      double cc = 0;
      double gs = 0;
      double gc = 0;
      for (int i=0;i<n;i++){
        sk->at(i) = (s->at(i) + kkk*z->at(i)*c->at(i));
        ck->at(i) = (c->at(i) - kkk*z->at(i)*s->at(i));
        ss += sk->at(i)*sk->at(i);
        sc += sk->at(i)*ck->at(i);
        cc += ck->at(i)*ck->at(i);
        gs += g->at(i)*sk->at(i);
        gc += g->at(i)*ck->at(i);
      }
      const double det = ss * cc - sc * sc;
      if (std::abs(det) < 1e-12) {
        return 1e12;
      }
      tx = (gc * sc - gs * cc) / det;
      ty = (gc * ss - gs * sc) / det;
      double sum = 0;
      for (int i=0;i<n;i++){
        sum+=(tx*sk->at(i)-ty*ck->at(i)+g->at(i))*(ty*z->at(i)*s->at(i)+tx*z->at(i)*c->at(i));
      }
      return sum;
    }
  };

  ROOT::Math::RootFinder& parabolicRF() {
    static thread_local ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kGSL_BISECTION);
    return rf;
  }

  struct IslCmp {
    bool operator() (const ActsExamples::IndexSourceLink& lhs,
                     const ActsExamples::IndexSourceLink& rhs) const {
      if (lhs.geometryId() < rhs.geometryId()) return true;
      if (rhs.geometryId() < lhs.geometryId()) return false;
      return lhs.index() < rhs.index();
    }
  };
}

ActsExamples::MySpacePointMaker::MySpacePointMaker(Config cfg, Acts::Logging::Level lvl)
  : IAlgorithm("MySpacePointMaker", lvl), m_cfg(std::move(cfg)) {
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);

  auto spConstructor =
      [](const Acts::Vector3& pos, std::optional<double> t,
         const Acts::Vector2& cov, std::optional<double> varT,
         boost::container::static_vector<Acts::SourceLink, 2> slinks)
      -> SimSpacePoint {
    return SimSpacePoint(pos, t, cov[0], cov[1], varT, std::move(slinks));
  };

  auto spBuilderConfig = Acts::SpacePointBuilderConfig();
  spBuilderConfig.trackingGeometry = m_cfg.detector->GetTrackingGeometry();
  m_slSurfaceAccessor.emplace(IndexSourceLink::SurfaceAccessor{*spBuilderConfig.trackingGeometry});
  spBuilderConfig.slSurfaceAccessor.connect<&IndexSourceLink::SurfaceAccessor::operator()>(&m_slSurfaceAccessor.value());
  m_spacePointBuilder = Acts::SpacePointBuilder<SimSpacePoint>(spBuilderConfig, spConstructor, Acts::getDefaultLogger("SpacePointBuilder", lvl));
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "misc-no-recursion"
ActsExamples::ProcessCode ActsExamples::MySpacePointMaker::execute(const AlgorithmContext& ctx) const {
  ACTS_INFO("Starting my space point maker for event " << ctx.eventNumber);
  using ISL = IndexSourceLink;
  using FtdLayerTypes = MyFtdGeo::FtdLayerTypes;

  std::shared_ptr<MyFtdGeo> ftdGeo = m_cfg.detector->FtdGeo();
  std::shared_ptr<MyFtdDetector> det = m_cfg.detector;

  const auto& measurements = m_inputMeasurements(ctx);
  const auto& measurementParticlesMap = m_inputMeasurementParticlesMap(ctx);

  std::vector<int> particleIds(measurements.size());
  for (auto meas : measurements) {
    for (const auto& [_, barcode] : ActsExamples::makeRange(measurementParticlesMap.equal_range(meas.index()))) {
      particleIds[meas.index()] = barcode.particle();
    }
  }

  // function to access measurement parameters using source links
  auto accessor = [&measurements](Acts::SourceLink slink) {
    const auto islink = slink.get<ISL>();
    const ConstVariableBoundMeasurementProxy meas = measurements.getMeasurement(islink.index());
    return std::make_pair(meas.fullParameters(), meas.fullCovariance());
  };

  // fill front and back strip vectors of pairs(source link, pair<stripEnd1, stripEnd2>)
  std::vector<std::pair<Acts::SourceLink, std::pair<Acts::Vector3, Acts::Vector3>>> frontStrips;
  std::vector<std::pair<Acts::SourceLink, std::pair<Acts::Vector3, Acts::Vector3>>> backStrips;
  std::vector<Acts::SourceLink> twoDimMeasurements;

  int nStations = ftdGeo->GetNumberOfStations();
  int nLayPerSt = ftdGeo->GetNLayersPerStation();
  int nLayers = ftdGeo->GetNumberOfLayers();

  std::map<ISL, std::tuple<Acts::Vector3, double, double, double>, IslCmp> islToXyzRot;

  std::vector<std::list<ISL>> mesPerLay(nLayers);
  for (const auto& isl : measurements.orderedIndices()) {
    const auto geoId = isl.geometryId();
    int iLayer = det->GeoIdToFtdLayer(geoId);
    mesPerLay[iLayer].emplace_back(isl);
    Acts::SourceLink slink{isl};
    // tube center shifted by the measured distance to wire
    const Acts::Surface* surface = m_slSurfaceAccessor.value()(slink);
    const auto [par, cov] = accessor(slink);
    auto xyz = surface->localToGlobal(ctx.geoContext, Acts::Vector2(par[0],0), Acts::Vector3());
    auto rot = surface->transform(ctx.geoContext).rotation();
    islToXyzRot[isl] = std::make_tuple(xyz, rot(0,1), rot(0,0), cov(0,0));
  }

  for (int iL = 0; iL < nLayers; ++iL) {
    auto& mesList = mesPerLay[iL];
    mesList.sort([](const ISL& a, const ISL& b) { return a.geometryId() < b.geometryId(); });
  }

  // iterate through layer measurements recursively
  auto constructCands = [&](auto&& self,
                            int layerId,
                            int station,
                            int maxDeltaStrawId,
                            int minMeasPerCand,
                            std::list<Candidate>& cands,
                            Candidate& cand,
                            ISL const* prevIsl) -> void
  {
    // not enough measurements found and not enough layers left to check
    if ((cand.sourceLinks.size() < minMeasPerCand) &&
        (nLayPerSt - (layerId % nLayPerSt) < (minMeasPerCand - static_cast<int>(cand.sourceLinks.size())))) {
      return;
    }

    // reached end of layers or another station -> store completed candidate
    if (ftdGeo->GetLayerStation(layerId) != station || layerId >= nLayers) {
      if (cand.sourceLinks.size() >= minMeasPerCand)
        cands.emplace_back(cand);
      return;
    }

    // skip fake pixel layers
    if (ftdGeo->GetLayerType(layerId) == FtdLayerTypes::kPixel) {
      self(self, layerId + 1, station, maxDeltaStrawId, minMeasPerCand, cands, cand, prevIsl);
      return;
    }

    auto& curLayerList = mesPerLay[layerId];

    if (prevIsl == nullptr) {
      // first layer we touch for this candidate -> try all measurements on this layer
      for (ISL& isl : curLayerList) {
        cand.sourceLinks.push_back(isl);
        self(self, layerId + 1, station, maxDeltaStrawId, minMeasPerCand, cands, cand, &cand.sourceLinks.back());
        cand.sourceLinks.pop_back(); // go back
      }
    } else {
      // iterate measurements in vicinity: [centerStrawId - maxDeltaStrawId, centerStrawId + maxDeltaStrawId]
      auto centerStrawId = prevIsl->geometryId().sensitive();
      auto it = curLayerList.begin();
      while (it != curLayerList.end() && it->geometryId().sensitive() < centerStrawId - maxDeltaStrawId) {
        ++it;
      }
      for (; it != curLayerList.end(); ++it) {
        auto strawId = it->geometryId().sensitive();
        if (strawId > centerStrawId + maxDeltaStrawId)
          break;
        cand.sourceLinks.push_back(*it);
        self(self, layerId + 1, station, maxDeltaStrawId, minMeasPerCand, cands, cand, &cand.sourceLinks.back());
        cand.sourceLinks.pop_back(); // go back
      }
      self(self, layerId + 1, station, maxDeltaStrawId, minMeasPerCand, cands, cand, prevIsl);
    }
  };

  // construct straw spacepoints
  std::list<Candidate> candidates[nStations];
  for (int iL = 0; iL < nLayers; ++iL) {
    int layerType = ftdGeo->GetLayerType(iL);
    if (layerType == FtdLayerTypes::kPixel) continue;
    int station = ftdGeo->GetLayerStation(iL);
    // if (station > 0) continue;
    if (station==1 || station==3) continue;
    auto& candList = candidates[station];
    Candidate cand;  // starts empty
    cand.station = station;
    int maxDeltaStrawId;
    if (station == 0) maxDeltaStrawId = m_cfg.maxDeltaStrawId1;
    if (station == 2) maxDeltaStrawId = m_cfg.maxDeltaStrawId2;
    if (station == 4) maxDeltaStrawId = m_cfg.maxDeltaStrawId3;
    constructCands(constructCands, iL, station, maxDeltaStrawId, m_cfg.minMeasPerCand, candList, cand, nullptr);
  }

  // construct everything else
  double rMax = ftdGeo->GetLayerRMax(ftdGeo->GetNumberOfLayers()-1);
  for (auto& isl : measurements.orderedIndices()) {
    const auto geoId = isl.geometryId();
    const auto volumeId = geoId.volume();
    const auto layerId = geoId.layer();
    int iLayer = det->GeoIdToFtdLayer(geoId);
    int iStation = ftdGeo->GetLayerStation(iLayer);
    Acts::SourceLink slink{isl};
    int layerType = ftdGeo->GetLayerType(iLayer);
    if (layerType == FtdLayerTypes::kPixel) {
      twoDimMeasurements.emplace_back(slink);
    } else if (0) {
      const auto [par, cov] = accessor(slink);
      const Acts::Surface *surface = m_slSurfaceAccessor.value()(slink);
      // TODO: more realistic strip dimensions including inner radii and half-station splitting
      auto gpos1 = surface->localToGlobal(ctx.geoContext, Acts::Vector2(par[0], -rMax), Acts::Vector3());
      auto gpos2 = surface->localToGlobal(ctx.geoContext, Acts::Vector2(par[0], rMax), Acts::Vector3());
      if (layerId % 4 == 2) {
        frontStrips.emplace_back(slink, std::make_pair(gpos1, gpos2));
      } else if (layerId % 4 == 0) {
        backStrips.emplace_back(slink, std::make_pair(gpos1, gpos2));
      }
    }
  }

  // return ActsExamples::ProcessCode::SUCCESS;
  ACTS_VERBOSE("making space points from straws");

  SimSpacePointContainer spacePoints;

  std::vector<double> c;
  std::vector<double> s;
  std::vector<double> z;
  std::vector<double> g;
  std::vector<double> d;

  for (int iStation=0;iStation<nStations;iStation++) {
    ACTS_VERBOSE("Station " << iStation << ": number of filtered candidates " << candidates[iStation].size());

    for (auto it = candidates[iStation].begin(); it != candidates[iStation].end();) {
      auto& candidate = *it;
      ACTS_VERBOSE("  candidate.size=" << candidate.sourceLinks.size());
      int n = candidate.sourceLinks.size();
      c.resize(n);
      s.resize(n);
      z.resize(n);
      g.resize(n);
      d.resize(n);
      for (size_t i=0; i<n; i++){
        auto& isl = candidate.sourceLinks[i];
        const auto geoId = isl.geometryId();
        const auto layerId = geoId.layer();
        const auto strawId = geoId.sensitive();
        auto itXR = islToXyzRot.find(isl);
        if (itXR == islToXyzRot.end())
          continue;
        auto& cachedXyzParCov = itXR->second;
        auto xyz = std::get<0>(cachedXyzParCov);
        double x = xyz[0];
        double y = xyz[1];
        double cosp = std::get<1>(cachedXyzParCov);
        double sinp = std::get<2>(cachedXyzParCov);
        z[i] = xyz[2];
        s[i] = z[i]*sinp;
        c[i] = z[i]*cosp;
        g[i] = -x*sinp + y*cosp;
        d[i] = sqrt(std::get<3>(cachedXyzParCov));
      }
      double txl, tyl, varxx, varyy, varxy;
      double chi2lin = linear(s, c, g, d, txl, tyl, varxx, varyy, varxy, 0);
      ACTS_VERBOSE("lin:  n=" << n << " tx=" << txl <<" ty=" << tyl << " chi2/ndf=" << chi2lin/(n-2));
      // setting variances obtained from the linear fit
      candidate.varxx = varxx;
      candidate.varxy = varxy;
      candidate.varyy = varyy;
      double tx, ty, k;
      double chi2par = parabolic(z, s, c, g, d, tx, ty, k, 0);
      double chi2ndf = n>3 ? chi2par/(n-3) : chi2par;
      ACTS_VERBOSE("par:  n=" << n << " tx=" << tx <<" ty=" << ty << " chi2/ndf=" << chi2ndf);
      // printf("chi2ndf=%f\n",chi2ndf);
      candidate.tx = tx;
      candidate.ty = ty;
      candidate.k = k;
      candidate.chi2 = chi2par;
      candidate.chi2ndf = chi2ndf;
      if (chi2ndf > m_cfg.maxChi2) {
        ACTS_VERBOSE("  erasing...");        
        it = candidates[iStation].erase(it);
      } else {
        ++it;
      }
    }

    ACTS_VERBOSE("Start filtering...");        
    ACTS_VERBOSE("Collecting number of candidates per measurement...");        
    int nCandidates = candidates[iStation].size();
    std::map<int,std::set<int>> candidatesPerStraw;
    std::vector<Candidate> selectedCandidates;    
    std::set<int> selectedCandidateIds;
    int i = 0;
    for (auto& candidate : candidates[iStation]){
      for (auto& isl : candidate.sourceLinks){
        int straw = isl.geometryId().layer()*10000+isl.geometryId().sensitive();
        candidate.straws.push_back(straw);
        candidatesPerStraw[straw].insert(i);
      }
      selectedCandidates.push_back(candidate);
      selectedCandidateIds.insert(i);
      i++;
    }

    ACTS_VERBOSE("Info on shared measurements...");
    for (int i=0;i<selectedCandidates.size();i++){
      auto& candidate = selectedCandidates[i];
      
      // printf("station %d: candidate %03d: chi2/ndf=%e particles: ", candidate.station, i, candidate.chi2ndf);
      // for (auto straw : candidate.straws) printf("%06d ", straw);
      // printf("\n");

      for (auto straw : candidate.straws){
        if (candidatesPerStraw[straw].size() > 1) {
          candidate.sharedStraws++;
        }
      }
      // printf("Candidate %d, meas: %d, shared: %d, chi2=%f\n", i, candidate.sourceLinks.size(), candidate.sharedStraws, candidate.chi2);
    }

    auto sharedStrawsComperator = [&selectedCandidates](std::size_t a, std::size_t b) {
      return selectedCandidates[a].sharedStraws < selectedCandidates[b].sharedStraws;
    };

    auto candidateComperator =  [&selectedCandidates](std::size_t a, std::size_t b) {
      int nMeasA = selectedCandidates[a].sourceLinks.size();
      int nMeasB = selectedCandidates[b].sourceLinks.size();
      auto relSharedA = 1.*selectedCandidates[a].sharedStraws/nMeasA;
      auto relSharedB = 1.*selectedCandidates[b].sharedStraws/nMeasB;
      if (relSharedA != relSharedB) return relSharedA < relSharedB;
      if (nMeasA == nMeasB) return selectedCandidates[a].chi2 < selectedCandidates[b].chi2;
      return nMeasA > nMeasB;
    };

    for (int i = 0; i < m_cfg.maximumIterations; i++) {
      if (selectedCandidateIds.empty()) break;
      auto maximumSharedStraws = *std::max_element(selectedCandidateIds.begin(), selectedCandidateIds.end(), sharedStrawsComperator);
      // printf("candidate %d with maximumSharedStraws=%d\n",maximumSharedStraws, selectedCandidates[maximumSharedStraws].sharedStraws);
      if (selectedCandidates[maximumSharedStraws].sharedStraws <= m_cfg.maximumSharedStraws) break;
      auto badCandidate = *std::max_element(selectedCandidateIds.begin(), selectedCandidateIds.end(), candidateComperator);
      // printf("badCandidate=%d %zu\n",badCandidate, selectedCandidates[badCandidate].sourceLinks.size());
      // remove bad track
      for (auto straw : selectedCandidates[badCandidate].straws) {
        candidatesPerStraw[straw].erase(badCandidate);
        if (candidatesPerStraw[straw].size()>1) continue;
        // only one candidate left
        auto it = *candidatesPerStraw[straw].begin();
        selectedCandidates[it].sharedStraws--;
      }
      selectedCandidateIds.erase(badCandidate);
    }
    
    ACTS_VERBOSE("Selected candidates...");
    for (auto& id : selectedCandidateIds) {
      auto& sp = selectedCandidates[id];
      if (sp.station==1 || sp.station==3) continue;
      int nlinks = sp.sourceLinks.size();
      printf("station %d: chi2/ndf=%e particles: ", sp.station, sp.chi2ndf);
      for (auto& isl : sp.sourceLinks) printf("%d ", particleIds[isl.index()]);
      printf("\n");
      Acts::SourceLink slink1{sp.sourceLinks[0]};
      Acts::SourceLink slink2{sp.sourceLinks[nlinks-1]};
      const Acts::Surface* surface = m_slSurfaceAccessor.value()(slink1);
      sp.z = surface->center(ctx.geoContext)[2];
      sp.x = sp.tx * sp.z + sp.k * sp.ty * sp.z * sp.z;
      sp.y = sp.ty * sp.z - sp.k * sp.tx * sp.z * sp.z;
      // printf("Candidate %d, meas: %d, shared: %d, chi2=%f\n", id, sp.sourceLinks.size(), sp.sharedStraws, sp.chi2);
      ACTS_VERBOSE("SP: k=" << sp.k << " x=" << sp.x << " y=" << sp.y << " z=" << sp.z << " var_xy=" << sp.varxy*sp.z*sp.z);

      boost::container::static_vector<Acts::SourceLink, 2> slinks = {slink1, slink2};
      Acts::Vector3 pos{sp.x,sp.y,sp.z};
      double var_r = (sp.x*sp.x*sp.varxx + sp.y*sp.y*sp.varyy + 2*sp.x*sp.y*sp.varxy)/(sp.x*sp.x+sp.y*sp.y);
      double var_z = 0.1;
      if (var_r>1e-5) continue;
      //std::back_inserter(spacePoints) = SimSpacePoint(pos, sp.chi2ndf*Acts::UnitConstants::ns, var_r, var_z, 0, slinks);
      std::back_inserter(spacePoints) = SimSpacePoint(pos, 0, var_r, var_z, 0, slinks);
      // std::back_inserter(spacePoints) = SimSpacePoint(pos, sp.varxy*sp.z*sp.z*Acts::UnitConstants::ns, sp.varxx*sp.z*sp.z, sp.varyy*sp.z*sp.z, 0, slinks);
    }
  }



  ACTS_VERBOSE("making strip pairs:" << "back: " << frontStrips.size() << " front: " << backStrips.size());

  // make space points from strips
  for (auto& fstrip : frontStrips) {
    float fz = fstrip.second.second[2];
    for (auto& bstrip : backStrips) {
      float bz = bstrip.second.second[2];
      if (fabs(fz - bz) > 20) continue;
      std::vector<Acts::SourceLink> slinks = {fstrip.first, bstrip.first};
      auto strippair = std::make_pair(fstrip.second, bstrip.second);
      Acts::SpacePointBuilderOptions spOptStrips{strippair, accessor};
      m_spacePointBuilder.buildSpacePoint(ctx.geoContext, slinks, spOptStrips, std::back_inserter(spacePoints));
    }
  }

  // make space points from 2D measurements
  Acts::SpacePointBuilderOptions spOpt;
  spOpt.paramCovAccessor = accessor;

  for (auto& slink : twoDimMeasurements) {
    // m_spacePointBuilder.buildSpacePoint(ctx.geoContext, {slink}, spOpt, std::back_inserter(spacePoints));
  }

  spacePoints.shrink_to_fit();

  ACTS_VERBOSE("Created " << spacePoints.size() << " space points");
  m_outputSpacePoints(ctx, std::move(spacePoints));

  return ActsExamples::ProcessCode::SUCCESS;
}

double ActsExamples::MySpacePointMaker::linear(
  std::vector<double> &s, std::vector<double> &c, std::vector<double> &g, std::vector<double> &d, 
  double &tx, double &ty, double &varxx, double &varyy, double &varxy, bool debug) const
{
  int n = s.size();
  double ss = 0.0;
  double cc = 0.0;
  double sc = 0.0;
  double gs = 0.0;
  double gc = 0.0;

  for (int i = 0; i < n; ++i){
    double w = 1.0 / (d[i] * d[i]);
    ss += w * s[i] * s[i];
    cc += w * c[i] * c[i];
    sc += w * s[i] * c[i];
    gs += w * g[i] * s[i];
    gc += w * g[i] * c[i];
  }
  double dd = cc * ss - sc * sc;
  tx = (sc * gc - cc * gs) / dd;
  ty = (ss * gc - sc * gs) / dd;
  varxx = cc / dd;
  varyy = ss / dd;
  varxy = sc / dd;

  double chi2 = 0;
  for (int i=0;i<n;i++){
    double chi = (s[i]*tx - c[i]*ty + g[i])/d[i];
    chi2 += chi*chi;
  }
  if (debug) printf("tx=%f ty=%f chi2=%f varxx=%e varyy=%e varxy=%e\n",tx, ty, chi2, varxx, varyy, varxy);

  return chi2;
}

double ActsExamples::MySpacePointMaker::parabolic(
  std::vector<double> &z, std::vector<double> &s, std::vector<double> &c, std::vector<double> &g,
  std::vector<double> &d, double &tx, double &ty, double &k, bool debug) const
{
  int n = s.size();
  std::vector<double> sk(n, 0.);
  std::vector<double> ck(n, 0.);

  FunctorFk fk;
  fk.n = n;
  fk.z = &z;
  fk.s = &s;
  fk.c = &c;
  fk.g = &g;
  fk.sk = &sk;
  fk.ck = &ck;

  double f0 = fk(0.);
  double dk = 1e-5;
  double dfdk = (fk(dk)-f0)/dk;
  double kk = -2*f0/dfdk;
  double fkk = fk(kk);
  if (!(f0 * fkk < 0.)) {
    return 1000.;
  }
  double kmin = kk>0 ? 0  : kk;
  double kmax = kk>0 ? kk :  0;
  if (debug) printf("%f %f %f %f\n",kmin, kmax, fk(kmin), fk(kmax));

  const int oldLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kBreak; // suppress warnings

  ROOT::Math::Functor1D functor(fk);
  ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kGSL_BISECTION);
  rf.SetFunction(functor, kmin, kmax);
  if (rf.Solve()) {
    k = rf.Root();
    double froot = fk(k);
    tx = fk.tx;
    ty = fk.ty;
    if (debug) printf("k=%e %f iterations=%d\n", k, fk(k), rf.Iterations());
  } else {
    return 1000.;
  }

  gErrorIgnoreLevel = oldLevel;

  if (debug) { 
    printf("checking tx=%f ty=%f k=%f\n",tx, ty, k);
    double sum_tx = 0;
    double sum_ty = 0;
    for (int i=0;i<n;i++){
      double a = s[i] + k*z[i]*c[i];
      double b = c[i] - k*z[i]*s[i];
      double d = a*tx - b*ty + g[i];
      sum_tx += d*a;
      sum_ty += d*b;
    }
    printf("sum_tx=%f\n",sum_tx);
    printf("sum_ty=%f\n",sum_ty);
  }

  double chi2 = 0;
  for (int i=0;i<n;i++){
    double chi = (sk[i]*tx - ck[i]*ty + g[i])/d[i];
    chi2 += chi*chi;
  }
  return chi2;
}
