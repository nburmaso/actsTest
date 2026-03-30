// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "MyDigitizationAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "ActsExamples/Digitization/ModuleClusters.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace ActsExamples {

MyDigitizationAlgorithm::MyDigitizationAlgorithm(Config config, Acts::Logging::Level level)
    : IAlgorithm("MyDigitizationAlgorithm", level), m_cfg(std::move(config)) {

  m_inputHits.initialize(m_cfg.inputSimHits);
  std::vector<std::pair<Acts::GeometryIdentifier, Digitizer>> digitizerInput;
  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
  m_outputMeasurementParticlesMap.initialize(m_cfg.outputMeasurementParticlesMap);
  m_outputMeasurementSimHitsMap.initialize(m_cfg.outputMeasurementSimHitsMap);
  m_outputParticleMeasurementsMap.initialize(m_cfg.outputParticleMeasurementsMap);
  m_outputSimHitMeasurementsMap.initialize(m_cfg.outputSimHitMeasurementsMap);

  for (std::size_t i = 0; i < m_cfg.digitizationConfigs.size(); ++i) {
    GeometricConfig geoCfg;
    Acts::GeometryIdentifier geoId = m_cfg.digitizationConfigs.idAt(i);
    const auto& digiCfg = m_cfg.digitizationConfigs.valueAt(i);
    switch (digiCfg.smearingDigiConfig.params.size()) {
      case 0u:
        digitizerInput.emplace_back(geoId, makeDigitizer<0u>(digiCfg));
        break;
      case 1u:
        digitizerInput.emplace_back(geoId, makeDigitizer<1u>(digiCfg));
        break;
      case 2u:
        digitizerInput.emplace_back(geoId, makeDigitizer<2u>(digiCfg));
        break;
      default:
        throw std::invalid_argument("Unsupported smearer size");
    }
  }

  m_digitizers = Acts::GeometryHierarchyMap<Digitizer>(digitizerInput);
}

ProcessCode MyDigitizationAlgorithm::execute(const AlgorithmContext& ctx) const {
  const auto& simHits = m_inputHits(ctx);
  // output containers
  MeasurementContainer measurements;
  IndexMultimap<SimBarcode> measurementParticlesMap;
  IndexMultimap<Index> measurementSimHitsMap;
  measurements.reserve(simHits.size());
  measurementParticlesMap.reserve(simHits.size());
  measurementSimHitsMap.reserve(simHits.size());
  // setup random number generator
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);

  for (const auto& simHitsGroup : groupByModule(simHits)) {
    Acts::GeometryIdentifier moduleGeoId = simHitsGroup.first;
    const auto& moduleSimHits = simHitsGroup.second;
    auto surfaceItr = m_cfg.surfaceByIdentifier.find(moduleGeoId);
    const Acts::Surface* surfacePtr = surfaceItr->second;
    auto digitizerItr = m_digitizers.find(moduleGeoId);
    std::visit([&](const auto& digitizer) {
      for (auto h = moduleSimHits.begin(); h != moduleSimHits.end(); ++h) {
        const auto& simHit = *h;
        const auto simHitIdx = simHits.index_of(h);
        auto boundParamsRes = Acts::transformFreeToBoundParameters(simHit.position(), simHit.time(), simHit.direction(), 0, *surfacePtr, ctx.geoContext, 1.e-4);
        if (!boundParamsRes.ok()) {
          ACTS_VERBOSE("!boundParamsRes.ok()");
          continue;
        }
        const auto& boundParams = *boundParamsRes;
        ACTS_VERBOSE("" << boundParams[0]);

        DigitizedParameters dParameters;
        auto res = digitizer.smearing(rng, simHit, *surfacePtr, ctx.geoContext);
        const auto& [par, cov] = res.value();
        for (Eigen::Index ip = 0; ip < par.rows(); ++ip) {
          dParameters.indices.push_back(digitizer.smearing.indices[ip]);
          dParameters.values.push_back(par[ip]);
          dParameters.variances.push_back(cov(ip, ip));
        }
        auto measurement = createMeasurement(measurements, moduleGeoId, dParameters);
        measurementParticlesMap.emplace_hint(measurementParticlesMap.end(), measurement.index(), simHits.nth(simHitIdx)->particleId());
        measurementSimHitsMap.emplace_hint(measurementSimHitsMap.end(), measurement.index(), simHitIdx);
        if (par.rows()==1 && m_cfg.doDuplicateStrawMeasurements) {
          DigitizedParameters dParameters2;
          for (Eigen::Index ip = 0; ip < par.rows(); ++ip) {
            dParameters2.indices.push_back(digitizer.smearing.indices[ip]);
            dParameters2.values.push_back(-par[ip]);
            dParameters2.variances.push_back(cov(ip, ip));
          }
          auto measurement = createMeasurement(measurements, moduleGeoId, dParameters2);
          measurementParticlesMap.emplace_hint(measurementParticlesMap.end(), measurement.index(), ActsFatras::Barcode());
          measurementSimHitsMap.emplace_hint(measurementSimHitsMap.end(), measurement.index(), simHitIdx);
        }
      }
    }, *digitizerItr);
  }

  ACTS_VERBOSE("measurements.size()=" << measurements.size());
  m_outputMeasurements(ctx, std::move(measurements));
  m_outputParticleMeasurementsMap(ctx, invertIndexMultimap(measurementParticlesMap));
  m_outputSimHitMeasurementsMap(ctx, invertIndexMultimap(measurementSimHitsMap));
  m_outputMeasurementParticlesMap(ctx, std::move(measurementParticlesMap));
  m_outputMeasurementSimHitsMap(ctx, std::move(measurementSimHitsMap));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
