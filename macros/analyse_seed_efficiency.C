#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "ActsFatras/EventData/Barcode.hpp"
#include "../MyFtdDetector.h"

bool isGoodFtd(int64_t layerMask, int minHits = 3){
  int nHits[5] = {0,0,0,0,0}; // number of hits per station
  for (int st=0;st<5;st++){
    for (int l=0;l<7;l++) {
      nHits[st] += ((layerMask & (1ull << (7*st+l))) > 0);
    }
  }
  if (nHits[0]<3) return 0;
  if (nHits[2]<3) return 0;
  if (nHits[4]<3) return 0;
  return 1;
}

void analyse_seed_efficiency(TString dir = "../build/test", double etaMean = 1.9, double etaDif = 0.05, bool seedable = 1){
//void analyse_seed_efficiency(TString dir = "../build/test", double etaMean = 1.9, double etaDif = 0.05, bool seedable = 1){
//void analyse_seed_efficiency(TString dir = "../build/notpc_pi_16_7deg", double etaMean = 1.6, double etaDif = 0.05, bool seedable = 1){
//void analyse_seed_efficiency(TString dir = "../acts_pi_16", double etaMean = 1.6, double etaDif = 0.05, bool seedable = 1){
//void analyse_seed_efficiency(TString dir = "../acts", double etaMean = 1.9, double etaDif = 0.05, bool seedable = 1){
//void analyse_seed_efficiency(TString dir = "../acts", double etaMean = 1.6, double etaDif = 0.05, bool seedable = 1){
  dir.Append("/");
  // int shift = 0; //  isroc = 0;  isframe = 0;
  int shift = 3; //  isroc = 1;  isframe = 1;

  // setup particles
  TFile* fPart = new TFile(TString(dir + "particles.root"));
  TTree* tPart = (TTree*) fPart->Get("particles");
  uint32_t part_event_id;
  auto part_id  = new std::vector<unsigned long>; 
  auto part_pdg = new vector<int>;
  auto part_vz  = new std::vector<float>;
  auto part_pt  = new std::vector<float>;
  auto part_eta = new std::vector<float>;
  tPart->SetBranchAddress("event_id",&part_event_id);
  tPart->SetBranchAddress("particle_type",&part_pdg);
  tPart->SetBranchAddress("vz",&part_vz);
  tPart->SetBranchAddress("pt",&part_pt);  
  tPart->SetBranchAddress("eta",&part_eta);  
  int nEvents = tPart->GetEntries();

  // setup seeds
  TFile* fSeeds = new TFile(TString(dir + "seeds.root"));
  TTree* tSeeds = (TTree*) fSeeds->Get("seeds");
  uint32_t seed_event_id;
  ULong64_t measId1;
  ULong64_t measId2;
  ULong64_t measId3;
  tSeeds->SetBranchAddress("event_id", &seed_event_id);
  tSeeds->SetBranchAddress("measurement_id_1", &measId1);
  tSeeds->SetBranchAddress("measurement_id_2", &measId2);
  tSeeds->SetBranchAddress("measurement_id_3", &measId3);
  
  // setup measurements
  TFile* fMeas = new TFile(TString(dir + "measurements.root"));
  TTree* tMeas = (TTree*) fMeas->Get("measurements");
  int32_t meas_event_id;
  int32_t meas_volume_id;
  int32_t meas_layer_id;
  vector<vector<uint32_t>> meas_particles; auto pmeas_particles = &meas_particles;
  tMeas->SetBranchAddress("event_nr",&meas_event_id);
  tMeas->SetBranchAddress("particles",&pmeas_particles);
  tMeas->SetBranchAddress("volume_id",&meas_volume_id);
  tMeas->SetBranchAddress("layer_id",&meas_layer_id);

  printf("map indices of first measurements per event + map fired layers per particle\n");
  std::vector<std::map<int,int64_t>> vEventParticleFtdLayerMask(nEvents);
  std::map<int,int> mapFirstMeaurementIndex;
  int previous_event_nr = -1;
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if (meas_event_id!=previous_event_nr) {
      if (meas_event_id%100==0) printf("Event=%d\n", meas_event_id);
      previous_event_nr = meas_event_id;
      mapFirstMeaurementIndex[meas_event_id] = im;
    }
    // if (meas_volume_id!=16) continue;  // with TPC volumes
    if (meas_volume_id!=6) continue;
    int layer = meas_layer_id - shift;
    auto& mapParticleFtdLayerMask = vEventParticleFtdLayerMask[meas_event_id];
    int ip = meas_particles[0][2] - 1;
    if (auto m = mapParticleFtdLayerMask.find(ip); m == mapParticleFtdLayerMask.end()) mapParticleFtdLayerMask[ip] = 0;
    mapParticleFtdLayerMask[ip] |= (1ull << layer);
  }

  TH1D* hMcPtPi = new TH1D("hMcPtPi","",100,0.,1.);
  TH1D* hMcPtPr = new TH1D("hMcPtPr","",100,0.,1.);
  TH1D* hRcPtPi = new TH1D("hRcPtPi","",100,0.,1.);
  TH1D* hRcPtPr = new TH1D("hRcPtPr","",100,0.,1.);

  printf("loop over mc particles\n");
  for (int ev=0;ev<tPart->GetEntries();ev++){
    tPart->GetEntry(ev);
    auto& mapParticleFtdLayerMask = vEventParticleFtdLayerMask[ev];
    for (int ip = 0; ip<part_pdg->size(); ip++){
      int64_t ftdLayerMask = mapParticleFtdLayerMask[ip];
      if (seedable && !isGoodFtd(ftdLayerMask)) continue;
      int pdg = part_pdg->at(ip);
      float vz = part_vz->at(ip);
      float pt = part_pt->at(ip);
      float eta = part_eta->at(ip);
      if (abs(eta-etaMean)<etaDif && abs(vz)<1.) {
        if (abs(pdg)== 211) hMcPtPi->Fill(pt);
        if (abs(pdg)==2212) hMcPtPr->Fill(pt);
      }
    }
  }

  printf("loop over seeds\n");  
  int previous_event_id = -1;
  for (int is=0;is<tSeeds->GetEntries();is++){
    tSeeds->GetEntry(is);
    if (previous_event_id!=seed_event_id){
      previous_event_id = seed_event_id;
      if (seed_event_id%100==0) printf("%d\n",seed_event_id);
      tPart->GetEntry(seed_event_id);
    }
    int mshift = mapFirstMeaurementIndex[seed_event_id];
    tMeas->GetEntry(mshift + measId1);
    auto hit_particle_id1 = meas_particles[0][2];
    tMeas->GetEntry(mshift + measId2);
    auto hit_particle_id2 = meas_particles[0][2];
    tMeas->GetEntry(mshift + measId3);
    auto hit_particle_id3 = meas_particles[0][2];
    // TODO loop over meas_particles (in case of several contributors to the measurement)
    if (hit_particle_id1!=hit_particle_id2 || hit_particle_id2!=hit_particle_id3) { // fake seed
      continue;
    }
    int ip = hit_particle_id1 - 1;
    int64_t ftdLayerMask = vEventParticleFtdLayerMask[seed_event_id][ip];
    if (seedable && !isGoodFtd(ftdLayerMask)) continue;    
    int pdg = part_pdg->at(ip);
    float vz = part_vz->at(ip);
    float pt = part_pt->at(ip);
    float eta = part_eta->at(ip);
    if (abs(eta-etaMean)<etaDif && abs(vz)<1.) {
      if (abs(pdg)== 211) hRcPtPi->Fill(pt);
      if (abs(pdg)==2212) hRcPtPr->Fill(pt);
    }
  }

  new TCanvas;
  hMcPtPi->Draw();

  new TCanvas;
  auto hEffPtPi = (TH1D*) hRcPtPi->Clone("hEffPtPi");
  hEffPtPi->Divide(hRcPtPi, hMcPtPi, 1, 1, "B");
  hEffPtPi->Draw();

  if (0) {
    new TCanvas;
    hMcPtPr->Draw();
 
    new TCanvas;
    auto hEffPtPr = (TH1D*) hRcPtPr->Clone("hEffPtPr");
    hEffPtPr->Divide(hRcPtPr, hMcPtPr, 1, 1, "B");
    hEffPtPr->Draw();
  }

  TFile* f = new TFile(TString(dir + "seed_efficiency.root"),"update");
  hMcPtPi->Write(Form("hMcPtPi%.0f",etaMean*10));
  hMcPtPr->Write(Form("hMcPtPr%.0f",etaMean*10));
  hRcPtPi->Write(Form("hRcPtPi%.0f",etaMean*10));
  hRcPtPr->Write(Form("hRcPtPr%.0f",etaMean*10));
  f->Close();

}
