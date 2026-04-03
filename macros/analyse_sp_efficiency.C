R__ADD_INCLUDE_PATH(/mnt/nvme0n1/nica/actsTest/build/stage/include)
R__ADD_LIBRARY_PATH(/mnt/nvme0n1/nica/actsTest/build/stage/lib)
R__LOAD_LIBRARY(libactsTestLib.so)
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "ActsFatras/EventData/Barcode.hpp"
#include "MyFtdGeo.h"
#include "MyFtdDetector.h"

const int nStations = 5;
const int nLayersPerStation = 9;
const int shift = 2;

int kPixel = 2;

std::shared_ptr<MyFtdDetector> det = nullptr;
std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
std::shared_ptr<MyFtdGeo> ftdGeo = nullptr;

bool isGoodSP(int64_t layerMask, int st){
  int nHits = 0; // number of hits per station
  for (int l=0;l<nLayersPerStation;l++) {
    if (ftdGeo->GetLayerType(l) == kPixel) continue;
    nHits += ((layerMask & (1ull << (nLayersPerStation*st+l))) > 0);
  }
  if (nHits<3) return 0;
  return 1;
}

void analyse_sp_efficiency(
  std::string inputDir = "./test/",
  double etaMean = 1.9, double etaDif = 0.1,
  int selected_station_sp = 0)
{
  det = std::make_shared<MyFtdDetector>();
  trackingGeometry = det->GetTrackingGeometry(true, true, true, false);
  ftdGeo = det->FtdGeo();

  // gStyle->SetOptStat(0);
  // setup particles
  TFile* fPart = new TFile(Form("%s/particles.root", inputDir.c_str()));
  TTree* tPart = (TTree*) fPart->Get("particles");
  UInt_t part_event_id;
  auto part_id  = new std::vector<uint32_t>; 
  auto part_pdg = new std::vector<int>;
  auto part_vz  = new std::vector<float>;
  auto part_px  = new std::vector<float>;
  auto part_py  = new std::vector<float>;
  auto part_pz  = new std::vector<float>;
  auto part_pt  = new std::vector<float>;
  auto part_eta = new std::vector<float>;
  auto part_phi = new std::vector<float>;
  auto part_mid = new std::vector<uint32_t>;
  tPart->SetBranchAddress("event_id",&part_event_id);
  tPart->SetBranchAddress("particle",&part_id);
  tPart->SetBranchAddress("particle_type",&part_pdg);
  tPart->SetBranchAddress("vz",&part_vz);
  tPart->SetBranchAddress("px",&part_px);
  tPart->SetBranchAddress("py",&part_py);
  tPart->SetBranchAddress("pz",&part_pz);  
  tPart->SetBranchAddress("pt",&part_pt);  
  tPart->SetBranchAddress("eta",&part_eta);  
  tPart->SetBranchAddress("phi",&part_phi);  
  tPart->SetBranchAddress("number_of_hits",&part_mid); // mother id

  int nEvents = tPart->GetEntries();

  // setup measurements
  TFile* fMeas = new TFile(Form("%s/measurements.root", inputDir.c_str()));
  TTree* tMeas = (TTree*) fMeas->Get("measurements");
  int32_t meas_event_id;
  int32_t meas_volume_id;
  int32_t meas_layer_id;
  vector<vector<uint32_t>> meas_particles; auto pmeas_particles = &meas_particles;
  tMeas->SetBranchAddress("event_nr",&meas_event_id);
  tMeas->SetBranchAddress("particles",&pmeas_particles);
  tMeas->SetBranchAddress("volume_id",&meas_volume_id);
  tMeas->SetBranchAddress("layer_id",&meas_layer_id);

  std::vector<std::vector<std::vector<uint32_t>>> vMeasParticleIds(nEvents);
  std::vector<std::map<int,int64_t>> vEventParticleFtdLayerMask(nEvents);
  printf("fill measurement particle ids + map fired layers per particle\n");
  int previous_event = -1;
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if (meas_event_id!=previous_event) {
      if (meas_event_id%100==0) printf("Event=%d\n", meas_event_id);
      previous_event = meas_event_id;
    }
    vMeasParticleIds[meas_event_id].push_back(vector<uint32_t>(meas_particles.size()));
    auto& mapParticleFtdLayerMask = vEventParticleFtdLayerMask[meas_event_id];
    for (int i=0;i<meas_particles.size();i++){
      int ip = meas_particles[i][2]-1;
      vMeasParticleIds[meas_event_id].back()[i] = ip;
      if (auto m = mapParticleFtdLayerMask.find(ip); m == mapParticleFtdLayerMask.end()) mapParticleFtdLayerMask[ip] = 0;
      mapParticleFtdLayerMask[ip] |= (1ull << (meas_layer_id-shift));
    }
  }

  TFile* fSpacepoints = new TFile(Form("%s/spacepoints.root", inputDir.c_str()));
  TTree* tSpacepoints = (TTree*) fSpacepoints->Get("spacepoints");
  tSpacepoints->Print();
  float sx;
  float sy;
  float sz;
  float st;  
  float varxx;
  float varxy;
  float varyy;
  UInt_t sevent_id;
  ULong64_t sgeometry_id;
  ULong64_t smeas_id, smeas_id_2;
  tSpacepoints->SetBranchAddress("x",&sx);
  tSpacepoints->SetBranchAddress("y",&sy);
  tSpacepoints->SetBranchAddress("z",&sz);
  tSpacepoints->SetBranchAddress("t",&varxy);
  tSpacepoints->SetBranchAddress("var_r",&varxx);  
  tSpacepoints->SetBranchAddress("var_z",&varyy);  
  tSpacepoints->SetBranchAddress("geometry_id",&sgeometry_id);
  tSpacepoints->SetBranchAddress("event_id",&sevent_id);
  tSpacepoints->SetBranchAddress("measurement_id", &smeas_id);
  tSpacepoints->SetBranchAddress("measurement_id_2", &smeas_id_2);

  std::map<int, std::map<int, int>> mSpoints;

  int previous_sevent_id = -1;
  for (int is=0;is<tSpacepoints->GetEntries();is++){
    tSpacepoints->GetEntry(is);
    if (previous_sevent_id!=sevent_id){
      previous_sevent_id = sevent_id;
      if (sevent_id%100==0) printf("%d\n",sevent_id);
    }
    int station = ftdGeo->GetLayerStation(det->GeoIdToFtdLayer(Acts::GeometryIdentifier(sgeometry_id)));
    if (station != selected_station_sp) continue;
    // printf("station %d\n", station);
    auto& partIds = vMeasParticleIds[sevent_id];
    for (auto id1 : partIds[smeas_id]) {
      for (auto id2 : partIds[smeas_id_2]) {
        if (id1!=id2) continue; // fake spoint
        mSpoints[sevent_id][id1] += 1;
      }
    }
  }

  TH1D* hSpRcPtPi   = new TH1D(Form("hSpRcPtPi%.0f_%d",   etaMean*10, selected_station_sp),"",100,0.,1.);
  TH1D* hSpRcPtPr   = new TH1D(Form("hSpRcPtPr%.0f_%d",   etaMean*10, selected_station_sp),"",100,0.,1.);
  TH1D* hSpablePtPi = new TH1D(Form("hSpablePtPi%.0f_%d", etaMean*10, selected_station_sp),"",100,0.,1.);
  TH1D* hSpablePtPr = new TH1D(Form("hSpablePtPr%.0f_%d", etaMean*10, selected_station_sp),"",100,0.,1.);

  for (int ev=0;ev<nEvents;ev++){ // events
    tPart->GetEntry(ev);
    for (int ip=0;ip<part_pdg->size();ip++){ // particles
      if (part_mid->at(ip)>0) continue; // only primaries
      auto& mapParticleFtdLayerMask = vEventParticleFtdLayerMask[ev];
      int64_t ftdLayerMask = mapParticleFtdLayerMask[ip];
      int pdg = part_pdg->at(ip);
      float vz = part_vz->at(ip);
      float pt = part_pt->at(ip);
      float eta = part_eta->at(ip);
      float phi = part_phi->at(ip);
      if (abs(eta-etaMean)>etaDif || abs(vz)>1.) continue;
      if (isGoodSP(ftdLayerMask, selected_station_sp)) {
        if (abs(pdg)== 211) hSpablePtPi->Fill(pt);
        if (abs(pdg)==2212) hSpablePtPr->Fill(pt);
        auto itMapSPointsEvent = mSpoints.find(ev);
        if (itMapSPointsEvent != mSpoints.end()) {
          auto& mapSPointsEvent = itMapSPointsEvent->second;
          auto itSpPart = mapSPointsEvent.find(ip);
          if (itSpPart != mapSPointsEvent.end()) {
            if (itSpPart->second > 0) {
              if (abs(pdg)== 211) hSpRcPtPi->Fill(pt);
              if (abs(pdg)==2212) hSpRcPtPr->Fill(pt);
            }
          }
        }
      }
    }
  }

  auto* c = new TCanvas;
  auto hEffSpPtPi = (TH1D*) hSpRcPtPi->Clone(Form("hEffSpPtPi%.0f_%d", etaMean*10., selected_station_sp));
  hEffSpPtPi->Divide(hSpRcPtPi, hSpablePtPi, 1, 1, "B");
  hEffSpPtPi->Draw();
  c->Print(Form("%s/sp_eff_pt_pi.png", inputDir.c_str()));

  auto* fout = new TFile(Form("%s/sp_efficiency.root", inputDir.c_str()), "update");
  hSpRcPtPi->Write(hSpRcPtPi->GetName(),TObject::kOverwrite);
  hSpRcPtPr->Write(hSpRcPtPr->GetName(),TObject::kOverwrite);
  hSpablePtPi->Write(hSpablePtPi->GetName(),TObject::kOverwrite);
  hSpablePtPr->Write(hSpablePtPr->GetName(),TObject::kOverwrite);
  fout->Close();
}
