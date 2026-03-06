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

ActsExamples::MySpacePointMaker::MySpacePointMaker(Config cfg, Acts::Logging::Level lvl)
  : IAlgorithm("MySpacePointMaker", lvl), m_cfg(std::move(cfg)) {
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);

  auto spConstructor =
      [](const Acts::Vector3& pos, std::optional<double> t,
         const Acts::Vector2& cov, std::optional<double> varT,
         boost::container::static_vector<Acts::SourceLink, 2> slinks)
      -> SimSpacePoint {
    return SimSpacePoint(pos, t, cov[0], cov[1], varT, std::move(slinks));
  };

  auto spBuilderConfig = Acts::SpacePointBuilderConfig();
  spBuilderConfig.trackingGeometry = m_cfg.trackingGeometry;
  m_slSurfaceAccessor.emplace(IndexSourceLink::SurfaceAccessor{*m_cfg.trackingGeometry});
  spBuilderConfig.slSurfaceAccessor.connect<&IndexSourceLink::SurfaceAccessor::operator()>(&m_slSurfaceAccessor.value());
  m_spacePointBuilder = Acts::SpacePointBuilder<SimSpacePoint>(spBuilderConfig, spConstructor, Acts::getDefaultLogger("SpacePointBuilder", lvl));
}

ActsExamples::ProcessCode ActsExamples::MySpacePointMaker::execute(const AlgorithmContext& ctx) const {
  const auto& measurements = m_inputMeasurements(ctx);

  // function to access measurement parameters using source links
  auto accessor = [&measurements](Acts::SourceLink slink) {
    const auto islink = slink.get<IndexSourceLink>();
    const ConstVariableBoundMeasurementProxy meas = measurements.getMeasurement(islink.index());
    return std::make_pair(meas.fullParameters(), meas.fullCovariance());
  };

  // fill front and back strip vectors of pairs(source link, pair<stripEnd1, stripEnd2>)
  std::vector<std::pair<Acts::SourceLink, std::pair<Acts::Vector3, Acts::Vector3>>> frontStrips;
  std::vector<std::pair<Acts::SourceLink, std::pair<Acts::Vector3, Acts::Vector3>>> backStrips;
  std::vector<Acts::SourceLink> twoDimMeasurements;

  int shift = 3; // FIXME
  int nStations = 5; // TODO read from config
//  std::list<std::vector<ActsExamples::IndexSourceLink>> candidates[nStations];
  std::list<Candidate> candidates[nStations];

  double rMax = 1300; // FIXME Extract from surface dimensions
  for (auto& isl : measurements.orderedIndices()) {
    const auto geoId = isl.geometryId();
    const auto volumeId = geoId.volume();
    const auto layerId = geoId.layer();
    int iStation = (layerId-shift)/7;
    int iLayer = (layerId-shift)%7;

    ACTS_DEBUG("volumeId= " << volumeId << " layerId=" << layerId);
    Acts::SourceLink slink{isl};

    if ((layerId-shift)%7==3 || layerId==2 || layerId==38) {
    // if (layerType[layerId - shift]==2) {
      twoDimMeasurements.emplace_back(slink);
    } else if (0){
      const auto [par, cov] = accessor(slink);
      const Acts::Surface* surface = m_slSurfaceAccessor.value()(slink);
      // TODO: more realistic strip dimensions including inner radii and half-station splitting
      auto gpos1 = surface->localToGlobal(ctx.geoContext, Acts::Vector2(par[0],-rMax), Acts::Vector3());
      auto gpos2 = surface->localToGlobal(ctx.geoContext, Acts::Vector2(par[0], rMax), Acts::Vector3());
      if (layerId%4==2) {
        frontStrips.emplace_back(std::make_pair(slink, std::make_pair(gpos1, gpos2)));
      } else if (layerId%4==0) {
        backStrips.emplace_back(std::make_pair(slink, std::make_pair(gpos1, gpos2)));
      }
    } else {
      int iStraw = geoId.sensitive();
      // construct sp candidates as vectors of index source links from the same station corresponding to
      // straws with similar straw ids (~similar angles)
      for (auto& candidate : candidates[iStation]){
        bool isCompatibleCandidate = 1;
        for (auto& isl2 : candidate.sourceLinks){
          const auto geoId2 = isl2.geometryId();
          int iLayer2 = (geoId2.layer()-shift)%7;
          int iStraw2 = geoId2.sensitive();
          int layerDif = iLayer - iLayer2;
          int strawDif = iStraw - iStraw2;
          if (layerDif<=0) { isCompatibleCandidate = 0; break ; } // should be imposible by construction
          if (layerDif==4) { // same direction
            if (strawDif<-1 || strawDif>0) { isCompatibleCandidate = 0; break ; }
          } else if (layerDif==1){
            if ((iLayer2==0 || iLayer2==4) && (strawDif<-2 || strawDif>3)) { isCompatibleCandidate = 0; break ; };
            if ((iLayer2==1 || iLayer2==5) && (strawDif<-6 || strawDif>4)) { isCompatibleCandidate = 0; break ; };
          } else if (layerDif==2){
            if ((iLayer2==0 || iLayer2==4) && (strawDif<-4 || strawDif>2)) { isCompatibleCandidate = 0; break ; };
          }
        }
        if (!isCompatibleCandidate) continue;
        candidate.sourceLinks.push_back(isl);
        ACTS_DEBUG("Adding new measurement to existing candidate");
      }
      if (iLayer==0 || iLayer==1 || iLayer==2 || iLayer==4){
        std::vector<ActsExamples::IndexSourceLink> new_candidate;
        new_candidate.push_back(isl);
        candidates[iStation].emplace_back(new_candidate,iStation);
        ACTS_DEBUG("Adding new candidate");
      }
    }
  }

  for (int iStation=0;iStation<nStations;iStation++) {
    ACTS_DEBUG("Station " << iStation << ": number of candidates " << candidates[iStation].size());    
    // for (auto candidate : candidates[iStation]) {
    for (auto it = candidates[iStation].begin(); it != candidates[iStation].end();) {
      auto& candidate = *it;
      ACTS_DEBUG("  candidate.size=" << candidate.sourceLinks.size());
      for (auto& isl : candidate.sourceLinks) {
        ACTS_DEBUG("   layer=" << isl.geometryId().layer() << " straw=" << isl.geometryId().sensitive());
      }
      if (candidate.sourceLinks.size()<3) {
        ACTS_DEBUG(  "erasing...");
        it = candidates[iStation].erase(it);
      } else {
        it++;
      }
    }
  }
  ACTS_DEBUG("making space points from straws");


  SimSpacePointContainer spacePoints;

  std::vector<double> c;
  std::vector<double> s;
  std::vector<double> z;
  std::vector<double> g;      
  std::vector<double> d;

  for (int iStation=0;iStation<nStations;iStation++) {
    ACTS_DEBUG("Station " << iStation << ": number of filtered candidates " << candidates[iStation].size());

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
        Acts::SourceLink slink{isl};
        const auto [par, cov] = accessor(slink);
        // tube center shifted by the measured distance to wire
        const Acts::Surface* surface = m_slSurfaceAccessor.value()(slink);
        auto xyz = surface->localToGlobal(ctx.geoContext, Acts::Vector2(par[0],0), Acts::Vector3());
        auto rot = surface->transform(ctx.geoContext).rotation();
        double x = xyz[0];
        double y = xyz[1];
        double cosp = rot(0,1);
        double sinp = rot(0,0);
        
        z[i] = xyz[2];
        s[i] = z[i]*sinp;
        c[i] = z[i]*cosp;
        g[i] = -x*sinp + y*cosp;
        d[i] = sqrt(cov(0,0));
        // ACTS_DEBUG("  x=" << x << " y=" << y <<" z=" << z[i]);  
      }
      double t,k,dt,dk;
      double chi2 = analytic(s, c, g, d, t, k, dt, dk);
      ACTS_DEBUG("  n=" << n << " t=" << t <<" k=" << k << " chi2/ndf=" << chi2/(n-2));
      double tx,ty,p;

      double chi2helix = helix(z, s, c, g, d, tx, ty, p);
//      ACTS_DEBUG("  n=" << n << " tx=" << tx <<" ty=" << ty << " chi2/ndf=" << chi2helix/(n-2));
      ACTS_DEBUG("  n=" << n << " t=" << sqrt(tx*tx+ty*ty) <<" k=" << ty/tx << " chi2/ndf=" << (n>3 ? chi2helix/(n-3) : chi2helix));
      candidate.tx = tx;
      candidate.ty = ty;
      candidate.k = p;
      // candidate.tx = t*cos(atan(k));
      // candidate.ty = t*sin(atan(k));
      // candidate.k = 0;

      candidate.chi2 = chi2helix;
//       double chi2parabolic = parabolic(z, s, c, g, d, tx, ty, p);
// //      ACTS_DEBUG("  n=" << n << " tx=" << tx <<" ty=" << ty << " chi2/ndf=" << chi2helix/(n-2));
//       ACTS_DEBUG("  n=" << n << " t=" << sqrt(tx*tx+ty*ty) <<" k=" << ty/tx << " chi2/ndf=" << (n>3 ? chi2parabolic/(n-3) : chi2parabolic));
      if (n>3) chi2helix/=(n-3);
      if (chi2helix>2000) {
        ACTS_DEBUG("  erasing...");        
        it = candidates[iStation].erase(it);
      } else {
        it++;
      }
    }

    ACTS_DEBUG("Start filtering...");        
    ACTS_DEBUG("Collecting number of candidates per measurement...");        
    int nCandidates = candidates[iStation].size();
    std::map<int,std::set<int>> candidatesPerMeasurement;
    std::vector<Candidate> selectedCandidates;    
    std::set<int> selectedCandidateIds;
    int i = 0;
    for (auto& candidate : candidates[iStation]){
      for (auto& isl : candidate.sourceLinks){
        int measId = isl.index();
        candidatesPerMeasurement[measId].insert(i);
      }
      selectedCandidates.push_back(candidate);
      selectedCandidateIds.insert(i);
      i++;
    }

    ACTS_DEBUG("Info on shared measurements...");        
    for (int i=0;i<selectedCandidates.size();i++){
      auto& candidate = selectedCandidates[i];
      for (auto& isl : candidate.sourceLinks){
        int measId = isl.index();
        if (candidatesPerMeasurement[measId].size() > 1) {
          candidate.sharedMeasurements++;
        }
      }
      printf("Candidate %d, meas: %d, shared: %d, chi2=%f\n", i, candidate.sourceLinks.size(), candidate.sharedMeasurements, candidate.chi2);
    }

    auto sharedMeasurementsComperator = [&selectedCandidates](std::size_t a, std::size_t b) {
      return selectedCandidates[a].sharedMeasurements < selectedCandidates[b].sharedMeasurements;
    };

    auto candidateComperator =  [&selectedCandidates](std::size_t a, std::size_t b) {
      int nMeasA = selectedCandidates[a].sourceLinks.size();
      int nMeasB = selectedCandidates[b].sourceLinks.size();
      auto relSharedA = 1.*selectedCandidates[a].sharedMeasurements/nMeasA;
      auto relSharedB = 1.*selectedCandidates[b].sharedMeasurements/nMeasB;
      if (relSharedA != relSharedB) return relSharedA < relSharedB;
      if (nMeasA == nMeasB) return selectedCandidates[a].chi2 < selectedCandidates[b].chi2;
      return nMeasA > nMeasB;
    };

    int maximumSharedHits = 1;
    for (int i = 0; i < 1000; i++) {
      if (selectedCandidateIds.empty()) break;
      auto maximumSharedMeasurements = *std::max_element(selectedCandidateIds.begin(), selectedCandidateIds.end(), sharedMeasurementsComperator);
      printf("candidate %d with maximumSharedMeasurements=%d\n",maximumSharedMeasurements, selectedCandidates[maximumSharedMeasurements].sharedMeasurements);
      if (selectedCandidates[maximumSharedMeasurements].sharedMeasurements <= maximumSharedHits) break;
      auto badCandidate = *std::max_element(selectedCandidateIds.begin(), selectedCandidateIds.end(), candidateComperator);
      printf("badCandidate=%d %zu\n",badCandidate, selectedCandidates[badCandidate].sourceLinks.size());
      // remove bad track
      for (auto& isl : selectedCandidates[badCandidate].sourceLinks) {
        int im = isl.index();
        candidatesPerMeasurement[im].erase(badCandidate);
        if (candidatesPerMeasurement[im].size()>1) continue;
        // only one candidate left
        auto it = *candidatesPerMeasurement[im].begin();
        selectedCandidates[it].sharedMeasurements--;
      }
      selectedCandidateIds.erase(badCandidate);
    }
    ACTS_DEBUG("Selected candidates...");        
    for (auto& id : selectedCandidateIds) {
      auto& sp = selectedCandidates[id];
      // TODO read from config
//      if (sp.station == 0) sp.z = 2070;
      if (sp.station == 0) sp.z = 2120;
      if (sp.station == 1) sp.z = 2340;
      if (sp.station == 2) sp.z = 2560;
      if (sp.station == 3) sp.z = 2780;
//      if (sp.station == 4) sp.z = 3050;
      if (sp.station == 4) sp.z = 3000;
      sp.x = sp.tx * sp.z + sp.k * sp.ty * sp.z * sp.z;
      sp.y = sp.ty * sp.z - sp.k * sp.tx * sp.z * sp.z;
      printf("Candidate %d, meas: %d, shared: %d, chi2=%f\n", id, sp.sourceLinks.size(), sp.sharedMeasurements, sp.chi2);
      ACTS_DEBUG("SP: x=" << sp.x << " y=" << sp.y << " z=" << sp.z);

      Acts::SourceLink slink1{sp.sourceLinks[0]};
      Acts::SourceLink slink2{sp.sourceLinks[sp.sourceLinks.size()-1]};
      //std::vector<Acts::SourceLink> slinks = {slink1, slink2};
      boost::container::static_vector<Acts::SourceLink, 2> slinks = {slink1, slink2};
      Acts::Vector3 pos{sp.x,sp.y,sp.z};
      std::back_inserter(spacePoints) = SimSpacePoint(pos, 0, 0, 0, 0, slinks);
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
    m_spacePointBuilder.buildSpacePoint(ctx.geoContext, {slink}, spOpt, std::back_inserter(spacePoints));
  }

  spacePoints.shrink_to_fit();

  ACTS_DEBUG("Created " << spacePoints.size() << " space points");
  m_outputSpacePoints(ctx, std::move(spacePoints));

  return ActsExamples::ProcessCode::SUCCESS;
}

/*
double ActsExamples::MySpacePointMaker::parabolic(
    const std::vector<double>& z,
    const std::vector<double>& s,
    const std::vector<double>& c,
    const std::vector<double>& g,
    const std::vector<double>& d,
    double &tx,
    double &ty,
    double &k
) const{
  const int n = z.size();
  Eigen::MatrixXd A(n, 4);
  Eigen::VectorXd b(n);
  for (int i = 0; i < n; ++i) {
    A(i, 0) =  s[i];           // t_x
    A(i, 1) = -c[i];           // t_y
    A(i, 2) =  z[i] * c[i];    // u_x
    A(i, 3) =  z[i] * s[i];    // u_y
    b(i) = -g[i];
  }

  // Solve least squares: min ||A*theta - b||
  Eigen::VectorXd theta = (A.transpose() * A).ldlt().solve(A.transpose() * b);

  tx = theta(0);
  ty = theta(1);
  double ux = theta(2);
  double uy = theta(3);

  double denom = tx * tx + ty * ty;
  if (denom == 0.0)
    throw std::runtime_error("Cannot extract curvature (zero slope)");

  k = (ux * tx + uy * ty) / denom;

  double chi2 = 0;
  for (int i=0;i<n;i++){
    double chi = (s[i]*tx - c[i]*ty + z[i]*c[i]*k*tx + z[i]*s[i]*k*ty + g[i])/d[i];
    chi2 += chi*chi;
  }
  return chi2;
}
*/

double ActsExamples::MySpacePointMaker::analytic(std::vector<double> &a, std::vector<double> &cz, std::vector<double> &g, std::vector<double> &s, 
                                                 double &t, double &k, double &dt, double &dk, bool debug) const{
  
  int n = a.size();
  std::vector<double> b(n,0);
  for (int i=0;i<n;i++) b[i] = -cz[i];
  std::vector<double> aa(n,0);
  std::vector<double> bb(n,0);
  double A = 0;
  double B = 0;
  double sigA2 = 0;
  double sigB2 = 0;
  double covAB = 0;

  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      double ab = (a[i]*b[j]-b[i]*a[j]);
      aa[i]+=b[j]*ab;
      bb[i]+=a[j]*ab;
    }
    A+= g[i]*aa[i];
    B+= g[i]*bb[i];
    sigA2+=s[i]*s[i]*aa[i]*aa[i];
    sigB2+=s[i]*s[i]*bb[i]*bb[i];
    covAB+=s[i]*s[i]*aa[i]*bb[i];
  }

  k = - B/A;
  dk = sqrt(B*B*sigA2+A*A*sigB2-2*A*B*covAB)/A/A;

  std::vector<double> abk(n,0);
  std::vector<double> vt(n,0);
  std::vector<double> dtdk(n,0);
  std::vector<double> dkdc(n,0);
  std::vector<double> dtdc(n,0);
  double num = 0;
  double den = 0;
  double p1 = 0;
  double p2 = 0;
  for (int i=0;i<n;i++){
    abk[i] = a[i]+b[i]*k;
    num += g[i]*abk[i];
    den += abk[i]*abk[i];
    p1 += a[i]*b[i];
    p2 += b[i]*b[i];
  }
  double cosa = 1./sqrt(1+k*k);
  double sina = k*cosa;

  t = -1./cosa*num/den;

  for (int i=0; i<n; i++){
    vt[i]  = -1./cosa*abk[i]/den;
    dtdk[i] = vt[i]*(k/(1+k*k) + b[i]/abk[i]-2*(p1+p2*k)/den);
    dkdc[i] = (aa[i]*B - bb[i]*A)/A/A;
  }

  double dt2 = 0;
  for (int i=0; i<n; i++){
    dtdc[i] = vt[i];
    for (int j=0; j<n; j++){
      dtdc[i] += g[j]*dtdk[j]*dkdc[i];
    }  

    dt2+=s[i]*s[i]*dtdc[i]*dtdc[i];
  }

  dt = sqrt(dt2);
  
  double chi2 = 0;
  for (int i=0;i<n;i++){
    double d = a[i]*t*cosa + b[i]*t*sina + g[i];
    chi2+=d*d/s[i]/s[i];
  }
  return chi2;
}

double ActsExamples::MySpacePointMaker::helix(
  std::vector<double> &z, std::vector<double> &s, std::vector<double> &c, std::vector<double> &g, 
  std::vector<double> &d, double &tx, double &ty, double &p, bool debug) const 
{
  int n = s.size();
  double sk[n]={0};
  double ck[n]={0};  
  auto fk = [&n,&s,&c,&z,&g,&tx,&ty,&sk,&ck](double k) {
    double ss = 0;
    double sc = 0;
    double cc = 0;
    double gs = 0;
    double gc = 0;
    for (int i=0;i<n;i++){
      sk[i] = (s[i] + k*z[i]*c[i]);
      ck[i] = (c[i] - k*z[i]*s[i]);
      ss += sk[i]*sk[i];
      sc += sk[i]*ck[i];
      cc += ck[i]*ck[i];
      gs += g[i]*sk[i];
      gc += g[i]*ck[i];
    }
    tx = (gc*sc - gs*cc)/(ss*cc - sc*sc);
    ty = (gc*ss - gs*sc)/(ss*cc - sc*sc);

    double sum = 0;
    for (int i=0;i<n;i++){
      sum+=(tx*sk[i]-ty*ck[i]+g[i])*(ty*z[i]*s[i]+tx*z[i]*c[i]);
    }
    return sum;
  };
  
  double f0 = fk(0);
  double dk = 1e-7;
  double dfdk = (fk(dk)-f0)/dk;
  double kk = -2*f0/dfdk;
  double kmin = kk>0 ? 0  : kk;
  double kmax = kk>0 ? kk :  0;
  // printf("%f %f %f %f\n",kmin, kmax, fk(kmin), fk(kmax));
  
  ROOT::Math::Functor1D functor(fk);
  ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kGSL_BISECTION);
  rf.SetFunction(functor, kmin, kmax);
  rf.Solve();
  p = rf.Root();
  printf("k=%e %f iterations=%d\n",p, fk(p), rf.Iterations());

  double sum_tx = 0;
  double sum_ty = 0;
  printf("check %f %f %f\n",tx, ty, p);
  for (int i=0;i<n;i++){
    double a = s[i] + p*z[i]*c[i];
    double b = c[i] - p*z[i]*s[i];
    double d = a*tx - b*ty + g[i];
    sum_tx += d*a;
    sum_ty += d*b;
  }
  printf("sum_tx=%f\n",sum_tx);
  printf("sum_ty=%f\n",sum_ty);

  double chi2 = 0;
  for (int i=0;i<n;i++){
    double chi = (sk[i]*tx - ck[i]*ty + g[i])/d[i];
    chi2 += chi*chi;
  }
  return chi2;

}
