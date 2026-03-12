#include "MyTrackWriter.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <TFile.h>
#include <TTree.h>

MyTrackWriter::MyTrackWriter(const MyTrackWriter::Config& config, Acts::Logging::Level level)
    : WriterT(config.inputTracks, "MyTrackWriter", level), m_cfg(config) {
  
  m_outputFile = new TFile(m_cfg.filePath.c_str(), "recreate");
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  m_outputTree->Branch("event_nr", &m_eventNr);
  m_outputTree->Branch("track_nr", &m_trackNr);
  m_outputTree->Branch("chi2Sum", &m_chi2Sum);
  m_outputTree->Branch("NDF", &m_NDF);
  m_outputTree->Branch("eLOC0_fit", &m_eLOC0_fit);
  m_outputTree->Branch("eLOC1_fit", &m_eLOC1_fit);
  m_outputTree->Branch("ePHI_fit", &m_ePHI_fit);
  m_outputTree->Branch("eTHETA_fit", &m_eTHETA_fit);
  m_outputTree->Branch("eQOP_fit", &m_eQOP_fit);
  m_outputTree->Branch("measurementIds", &m_measurementIds);
}

MyTrackWriter::~MyTrackWriter() {
  m_outputFile->Close();
}

ProcessCode MyTrackWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();
  return ProcessCode::SUCCESS;
}

ProcessCode MyTrackWriter::writeT(const AlgorithmContext& ctx, const ConstTrackContainer& tracks) {
  m_eventNr = ctx.eventNumber;

  for (const auto& track : tracks) {
    m_trackNr.push_back(track.index());
    m_chi2Sum.push_back(track.chi2());
    m_NDF.push_back(track.nDoF());
    const auto& param = track.parameters();
    m_eLOC0_fit.push_back(param[Acts::eBoundLoc0]);
    m_eLOC1_fit.push_back(param[Acts::eBoundLoc1]);
    m_ePHI_fit.push_back(param[Acts::eBoundPhi]);
    m_eTHETA_fit.push_back(param[Acts::eBoundTheta]);
    m_eQOP_fit.push_back(param[Acts::eBoundQOverP]);
    std::vector<int> measurementIds;
    for (const auto& state : track.trackStatesReversed()) {
      if (!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) continue;
      auto measId = state.getUncalibratedSourceLink().template get<ActsExamples::IndexSourceLink>().index();
      measurementIds.push_back(measId);
    }
    m_measurementIds.push_back(measurementIds);
  }
  m_outputTree->Fill();
  m_trackNr.clear();
  m_chi2Sum.clear();
  m_NDF.clear();
  m_eLOC0_fit.clear();
  m_eLOC1_fit.clear();
  m_ePHI_fit.clear();
  m_eTHETA_fit.clear();
  m_eQOP_fit.clear();
  m_measurementIds.clear();
  return ProcessCode::SUCCESS;
}
