R__ADD_INCLUDE_PATH(/home/ekryshen/mpd/actsTest/build/stage/include)
R__ADD_LIBRARY_PATH(/home/ekryshen/mpd/actsTest/build/stage/lib)
R__LOAD_LIBRARY(libactsTestLib.so)
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "ActsFatras/EventData/Barcode.hpp"
#include "MyFtdGeo.h"
#include "MyFtdDetector.h"

const int nStations = 5;
const int nLayersPerStation = 10;
const int shift = 2;
const int minMeasPerCand = 4;

int kPixel = 2;

std::shared_ptr<MyFtdDetector> det = nullptr;
std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
std::shared_ptr<MyFtdGeo> ftdGeo = nullptr;

bool isGoodSP(int64_t layerMask, int st){
  int nHits = 0; // number of hits per station
  int type4 = 0;
  int type5 = 0;
  int type6 = 0;
  for (int l=0;l<nLayersPerStation;l++) {
    int layerIndex = nLayersPerStation*st+l;
    int type = ftdGeo->GetLayerType(layerIndex);
    if (type == kPixel) continue;
    bool isHit = TESTBIT(layerMask, layerIndex);
    if (!isHit) continue;
    if (type == 4) type4++;
    if (type == 5) type5++;
    if (type == 6) type6++;
    nHits++;
  }
  if (type4 == 0) return 0;
  if (type5 == 0) return 0;
  if (type6 == 0) return 0;
  if (nHits<minMeasPerCand) return 0;
  return 1;
}

void analyse_sp_efficiency(
  std::string inputDir = "../build/test/",
  double etaMean = 1.75, double etaDif = 0.2,
  int selected_station_sp = 0)
{
  // gStyle->SetOptStat(0);

  #define axisPt 100,0.,1.
  #define axisPhi 90,-M_PI,M_PI
  #define axisEta 50,1.5,2.0

  TH1D* hSpMcPtPi   = new TH1D(Form("hSpMcPtPi%.0f_%d",   etaMean*10, selected_station_sp),"",axisPt);
  TH1D* hSpMcPtPr   = new TH1D(Form("hSpMcPtPr%.0f_%d",   etaMean*10, selected_station_sp),"",axisPt);
  TH1D* hSpRcPtPi   = new TH1D(Form("hSpRcPtPi%.0f_%d",   etaMean*10, selected_station_sp),"",axisPt);
  TH1D* hSpRcPtPr   = new TH1D(Form("hSpRcPtPr%.0f_%d",   etaMean*10, selected_station_sp),"",axisPt);
  TH1D* hSpablePtPi = new TH1D(Form("hSpablePtPi%.0f_%d", etaMean*10, selected_station_sp),"",axisPt);
  TH1D* hSpablePtPr = new TH1D(Form("hSpablePtPr%.0f_%d", etaMean*10, selected_station_sp),"",axisPt);
  
  TH1D* hSpMcEtaPi   = new TH1D(Form("hSpMcEtaPi%.0f_%d",   etaMean*10, selected_station_sp),"",axisEta);
  TH1D* hSpMcEtaPr   = new TH1D(Form("hSpMcEtaPr%.0f_%d",   etaMean*10, selected_station_sp),"",axisEta);
  TH1D* hSpableEtaPi = new TH1D(Form("hSpableEtaPi%.0f_%d", etaMean*10, selected_station_sp),"",axisEta);
  TH1D* hSpableEtaPr = new TH1D(Form("hSpableEtaPr%.0f_%d", etaMean*10, selected_station_sp),"",axisEta);

  det = std::make_shared<MyFtdDetector>();
  trackingGeometry = det->GetTrackingGeometry(true, true, true, false);
  ftdGeo = det->FtdGeo();

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
  std::vector<std::vector<int>> vSpoints(nEvents);
  std::vector<std::vector<int64_t>> vFtdLayerMask(nEvents);
  for (int ev=0;ev<nEvents;ev++){ // events
    tPart->GetEntry(ev);
    int nParts = part_pdg->size();
    vSpoints[part_event_id].resize(nParts+1,0);
    vFtdLayerMask[part_event_id].resize(nParts+1,0);
  }

  // setup measurements
  TFile* fMeas = new TFile(Form("%s/measurements.root", inputDir.c_str()));
  TTree* tMeas = (TTree*) fMeas->Get("measurements");
  int32_t meas_event_id;
  int32_t meas_volume_id;
  int32_t meas_layer_id;
  float meas_true_theta;
  vector<vector<uint32_t>> meas_particles; auto pmeas_particles = &meas_particles;
  tMeas->SetBranchAddress("event_nr",&meas_event_id);
  tMeas->SetBranchAddress("particles",&pmeas_particles);
  tMeas->SetBranchAddress("volume_id",&meas_volume_id);
  tMeas->SetBranchAddress("layer_id",&meas_layer_id);
  tMeas->SetBranchAddress("true_theta",&meas_true_theta);
  int previous_event = -1;

  TH1D* hLayers = new TH1D("hLayers","",100,0,100);
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if (meas_event_id!=previous_event) {
      if (meas_event_id%100==0) printf("Event=%d\n", meas_event_id);
      previous_event = meas_event_id;
    }
    for (int i=0;i<meas_particles.size();i++){
      int ip = meas_particles[i][2]; // counted from 1
      SETBIT(vFtdLayerMask[meas_event_id][ip], (meas_layer_id - shift));
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
  float smajority;
  UInt_t sevent_id;
  ULong64_t sgeometry_id;
  ULong64_t smeas_id, smeas_id_2;
  tSpacepoints->SetBranchAddress("x",&sx);
  tSpacepoints->SetBranchAddress("y",&sy);
  tSpacepoints->SetBranchAddress("z",&sz);
  tSpacepoints->SetBranchAddress("t",&smajority);
  tSpacepoints->SetBranchAddress("var_r",&varxx);  
  tSpacepoints->SetBranchAddress("var_z",&varyy);  
  tSpacepoints->SetBranchAddress("geometry_id",&sgeometry_id);
  tSpacepoints->SetBranchAddress("event_id",&sevent_id);
  tSpacepoints->SetBranchAddress("measurement_id", &smeas_id);
  tSpacepoints->SetBranchAddress("measurement_id_2", &smeas_id_2);

  int previous_sevent_id = -1;
  for (int is=0;is<tSpacepoints->GetEntries();is++){
    tSpacepoints->GetEntry(is);
    if (previous_sevent_id!=sevent_id){
      previous_sevent_id = sevent_id;
      if (sevent_id%100==0) printf("%d\n",sevent_id);
    }
    int station = ftdGeo->GetLayerStation(det->GeoIdToFtdLayer(Acts::GeometryIdentifier(sgeometry_id)));
    if (station != selected_station_sp) continue;
    int majorityId = trunc(smajority);
    float majFrac = (smajority - majorityId) * 10.;
    // printf("majorityId=%d majFrac=%f\n", majorityId, majFrac);
    if (majFrac > 0.74) {
      vSpoints[sevent_id][majorityId] += 1;
    }
  }

  for (int ev=0;ev<nEvents;ev++){ // events
    tPart->GetEntry(ev);
    for (int i=0;i<part_pdg->size();i++){ // particles
      int ip=i+1; // +1 since particleId is counted from 1
      if (part_mid->at(i)>0) continue; // only primaries
      int pdg = part_pdg->at(i);
      float vz = part_vz->at(i);
      float pt = part_pt->at(i);
      float eta = part_eta->at(i);
      float phi = part_phi->at(i);
      if (abs(eta-etaMean)>etaDif || abs(vz)>1.) continue;
      if (abs(pdg)== 211) hSpMcPtPi->Fill(pt);
      if (abs(pdg)==2212) hSpMcPtPr->Fill(pt);
      if (abs(pdg)== 211) hSpMcEtaPi->Fill(eta);
      if (abs(pdg)==2212) hSpMcEtaPr->Fill(eta);
      auto& ftdLayerMask = vFtdLayerMask[part_event_id][ip];
      if (!isGoodSP(ftdLayerMask, selected_station_sp)) continue;
      if (abs(pdg)== 211) hSpablePtPi->Fill(pt);
      if (abs(pdg)==2212) hSpablePtPr->Fill(pt);
      if (abs(pdg)== 211) hSpableEtaPi->Fill(eta);
      if (abs(pdg)==2212) hSpableEtaPr->Fill(eta);
      if (vSpoints[part_event_id][ip]==0) continue;
      if (abs(pdg)== 211) hSpRcPtPi->Fill(pt);
      if (abs(pdg)==2212) hSpRcPtPr->Fill(pt);
    }
  }

  float mc = hSpMcPtPi->Integral();
  float rc = hSpRcPtPi->Integral();
  float spable = hSpablePtPi->Integral();
  printf("spable/mc=%.0f/%.0f=%.4f\n",spable, mc, spable/mc);
  printf("rc/spable = %.0f/%.0f=%.4f\n",rc, spable, rc/spable);

  new TCanvas;
  auto hEffSpableEtaPi = (TH1D*) hSpableEtaPi->Clone(Form("hEffSpableEtaPi%.0f_%d", etaMean*10., selected_station_sp));
  hEffSpableEtaPi->Divide(hSpableEtaPi, hSpMcEtaPi, 1, 1, "B");
  hEffSpableEtaPi->Draw();

  new TCanvas;
  auto hEffSpablePtPi = (TH1D*) hSpablePtPi->Clone(Form("hEffSpablePtPi%.0f_%d", etaMean*10., selected_station_sp));
  hEffSpablePtPi->Divide(hSpablePtPi, hSpMcPtPi, 1, 1, "B");
  hEffSpablePtPi->Draw();
  
  new TCanvas;
  auto hEffSpPtPi = (TH1D*) hSpRcPtPi->Clone(Form("hEffSpPtPi%.0f_%d", etaMean*10., selected_station_sp));
  hEffSpPtPi->Divide(hSpRcPtPi, hSpablePtPi, 1, 1, "B");
  hEffSpPtPi->Draw();
  gPad->Print(Form("%s/sp_eff_pt_pi.png", inputDir.c_str()));

  auto* fout = new TFile(Form("%s/sp_efficiency.root", inputDir.c_str()), "update");
  hSpRcPtPi->Write(hSpRcPtPi->GetName(),TObject::kOverwrite);
  hSpRcPtPr->Write(hSpRcPtPr->GetName(),TObject::kOverwrite);
  hSpablePtPi->Write(hSpablePtPi->GetName(),TObject::kOverwrite);
  hSpablePtPr->Write(hSpablePtPr->GetName(),TObject::kOverwrite);
  hSpMcPtPi->Write(hSpMcPtPi->GetName(),TObject::kOverwrite);
  hSpMcPtPr->Write(hSpMcPtPr->GetName(),TObject::kOverwrite);
  fout->Close();
}
