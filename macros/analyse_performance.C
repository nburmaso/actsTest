#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "ActsFatras/EventData/Barcode.hpp"
#include "tree_summary.C"

bool isGoodFtd(int64_t layerMask, int minHits = 5){
  int nHits[5] = {0,0,0,0,0}; // number of hits per station
  for (int st=0;st<5;st++){
    for (int l=0;l<7;l++) {
      if (l==3) continue;
      nHits[st] += ((layerMask & (1ull << (7*st+l))) > 0);
    }
  }
  // printf("%d %d %d %d %d\n",nHits[0],nHits[1],nHits[2],nHits[3],nHits[4]);
  if (nHits[0]<3) return 0;
  if (nHits[1]<3) return 0;
  if (nHits[2]<3) return 0;
  if (nHits[3]<3) return 0;
  if (nHits[4]<3) return 0;
  return 1;
}

bool isGoodRecoFtd(int64_t layerMask, int minHits = 5){
  int nHits[5] = {0,0,0,0,0}; // number of hits per station
  for (int st=0;st<5;st++){
    for (int l=0;l<7;l++) {
      if (l==3) continue;
      nHits[st] += ((layerMask & (1ull << (7*st+l))) > 0);
    }
  }
  if (nHits[0]<1) return 0;
  if (nHits[1]<1) return 0;
  if (nHits[2]<1) return 0;
  if (nHits[3]<1) return 0;
  if (nHits[4]<1) return 0;
  return 1;
}

bool isGoodSeed(int64_t layerMask, int minHits = 5){
  int nSeeds[5] = {0,0,0,0,0}; // number of seeds per station
  for (int st=0;st<5;st++){
    for (int l=0;l<7;l++) {
      if (l==3) continue;
      nSeeds[st] += ((layerMask & (1ull << (7*st+l))) > 0);
    }
  }
  if (nSeeds[0]<3) return 0;
  if (nSeeds[2]<3) return 0;
  if (nSeeds[4]<3) return 0;
  return 1;
}


// bool isGoodSeed(int64_t layerMask, int minHits = 5, int shift = 3){
//   int nSeeds[5] = {0,0,0,0,0}; // number of seeds per station
//   int nSeedsB = (layerMask & (1ull << (shift -1))) > 0;
//   int nSeedsF = (layerMask & (1ull << (7*4 + shift + 7))) > 0;
//   for (int st=0;st<5;st++){
//     nSeeds[st] += (layerMask & (1ull << (7*st + shift + 3))) > 0;
//   }
//   // printf("%d %d %d\n", nSeedsB, nSeeds[2], nSeedsF);
//   if (nSeedsB<1) return 0;
//   if (nSeedsF<1) return 0;
//   // if (nSeeds[0]<1) return 0;
//   if (nSeeds[2]<1) return 0;
//   // if (nSeeds[4]<1) return 0;
//   return 1;
// }

void analyse_performance(TString dir = "../build/test/", double etaMean = 1.9, double etaDif = 0.1, bool refit = 0, bool trackable = 1){
//void analyse_performance(TString dir = "../acts/", double etaMean = 1.7, double etaDif = 0.1, bool refit = 0, bool trackable = 1){
//void analyse_performance(TString dir = "../pi90/acts/", double etaMean = 1.7, double etaDif = 0.1, bool refit = 0, bool trackable = 1){
//void analyse_performance(TString dir = "../acts_pi_16/", double etaMean = 1.6, double etaDif = 0.1, bool refit = 0, bool trackable = 1){
  // gStyle->SetOptStat(0);
  int shift = 3; //  isroc = 1;  isframe = 1;
  // setup particles
  TFile* fPart = new TFile(TString(dir + "particles.root"));
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

  printf("fill match array\n");
  TFile* fTrack = new TFile(TString(dir + (refit ? "trackrefit.root" : "tracksummary.root")) );
  TTree* tTrack = (TTree*) fTrack->Get("tracksummary");
  SetBranchAddresses(tTrack);
  vector<vector<int>> vMatched(nEvents);
  vector<vector<int>> vNumberOfMatched(nEvents);
  vector<vector<int>> vSeeds(nEvents);

  TH1D* hLayers = new TH1D("hLayers","",40,0,40);
  TH1D* hLOC0fit = new TH1D("hLOC0fit","loc0",100,-100,100);
  TH1D* hLOC1fit = new TH1D("hLOC1fit","loc1",100,-100,100);
  TH1D* hNLayers = new TH1D("hNLayers","layers",40,0,40);
  for (int ev=0; ev<nEvents; ev++){
    if (ev%100==0) printf("Event = %d\n",ev);
    tPart->GetEntry(ev);
    vNumberOfMatched[ev].assign(part_id->size(),0);
    vMatched[ev].assign(part_id->size(),0);
    vSeeds[ev].assign(part_id->size(),0);
    tTrack->GetEntry(ev);

    for (int it=0; it<m_majorityParticleId->size(); it++){
      auto barcode = ActsFatras::Barcode().withData(m_majorityParticleId->at(it));
      // auto loc0 = m_eLOC0_fit->at(it);
      // auto loc1 = m_eLOC1_fit->at(it);
      // hLOC0fit->Fill(loc0);
      // hLOC1fit->Fill(loc1);
      // if (fabs(loc0)>8) continue;
      // if (fabs(loc1)>8) continue;
      int ip = barcode.particle()-1;
      if (ip<0) continue;
      if (vMatched[ev][ip]<1) vMatched[ev][ip]=1;
      if (!m_hasFittedParams->at(it)) continue;
      if (vMatched[ev][ip]<2) vMatched[ev][ip]=2;
      auto& layers = m_measurementLayer->at(it);
      int64_t trackFtdLayerMask = 0;
      for (int il=0;il<layers.size();il++) trackFtdLayerMask |= (1ull << (layers[il]-shift));
      for (int il=0;il<layers.size();il++) hLayers->Fill(layers[il]);
      if (trackable && !isGoodRecoFtd(trackFtdLayerMask)) continue;
      if (vMatched[ev][ip]<3) vMatched[ev][ip]=3;
      vNumberOfMatched[ev][ip]++;
      hNLayers->Fill(layers.size());
    }
  }


  new TCanvas;
  hLayers->Draw();
  new TCanvas;
  hLOC0fit->Draw();
  new TCanvas;
  hLOC1fit->Draw();
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

  TH1D* hSeedPtPi = new TH1D("hSeedPtPi","",100,0.,1.);
  TH1D* hSeedPtPr = new TH1D("hSeedPtPr","",100,0.,1.);
  TH1D* hSeedPhiPi = new TH1D("hSeedPhiPi","",180,-M_PI,M_PI);
  TH1D* hSeedPhiPr = new TH1D("hSeedPhiPr","",180,-M_PI,M_PI);

  printf("loop over seeds\n");
  int previous_seed_event_id = -1;
  for (int is=0;is<tSeeds->GetEntries();is++){
    tSeeds->GetEntry(is);
    if (previous_seed_event_id!=seed_event_id){
      previous_seed_event_id = seed_event_id;
      if (seed_event_id%100==0) printf("%d\n",seed_event_id);
    }
    auto& partIds = vMeasParticleIds[seed_event_id];
    for (auto id1 : partIds[measId1]) {
      for (auto id2 : partIds[measId2]) {
        for (auto id3 : partIds[measId3]) {
          if (id1!=id2 || id2!=id3) continue; // fake seed
          vSeeds[seed_event_id][id1]+=1;
        }
      }
    }
  }

  TH1D* hRcPtPi = new TH1D("hRcPtPi","",100,0.,1.);
  TH1D* hRcPtPr = new TH1D("hRcPtPr","",100,0.,1.);
  TH1D* hRcPhiPi = new TH1D("hRcPhiPi","",180,-M_PI,M_PI);
  TH1D* hRcPhiPr = new TH1D("hRcPhiPr","",180,-M_PI,M_PI);
  TH1D* hNumberOfMatchedPi = new TH1D("hNumberOfMatchedPi","",100,0.,1.);
  TH1D* hNumberOfGoodPi = new TH1D("hNumberOfGoodPi","",100,0.,1.);

  // setup performance
  // TFile* fPerf = new TFile(TString(dir + "performance_ckf.root"));
  // TTree* tPerf = (TTree*) fPerf->Get("matchingdetails");
  // uint32_t perf_event_nr;
  // std::vector<uint32_t> perf_barcode; auto perf_particle_id = &perf_barcode;
  // bool perf_matched;
  // tPerf->Print();
  // tPerf->SetBranchAddress("event_nr",&perf_event_nr);
  // tPerf->SetBranchAddress("particle_id",&perf_particle_id);
  // tPerf->SetBranchAddress("matched",&perf_matched);

  for (int ev=0;ev<nEvents;ev++){ // events
    tPart->GetEntry(ev);
    for (int ip=0;ip<part_pdg->size();ip++){ // particles
      if (part_mid->at(ip)>0) continue; // only primaries
      auto& mapParticleFtdLayerMask = vEventParticleFtdLayerMask[ev];
      int64_t ftdLayerMask = mapParticleFtdLayerMask[ip];
      if (trackable && (!isGoodFtd(ftdLayerMask) || !isGoodSeed(ftdLayerMask))) continue;    
      int pdg = part_pdg->at(ip);
      float vz = part_vz->at(ip);
      float pt = part_pt->at(ip);
      float eta = part_eta->at(ip);
      float phi = part_phi->at(ip);    
      if (abs(eta-etaMean)>etaDif || abs(vz)>1.) continue;
      if (vSeeds[ev][ip]==0) continue;
      if (vSeeds[ev][ip]>1) printf("Warning: seeds = %d\n",vSeeds[ev][ip]);
      if (abs(pdg)== 211) hSeedPtPi->Fill(pt);
      if (abs(pdg)==2212) hSeedPtPr->Fill(pt);
      if (abs(pdg)== 211) hSeedPhiPi->Fill(phi);
      if (abs(pdg)==2212) hSeedPhiPr->Fill(phi);
      if (vMatched[ev][ip]<3) continue;
      if (abs(pdg)== 211) hRcPtPi->Fill(pt);
      if (abs(pdg)==2212) hRcPtPr->Fill(pt);
      if (abs(pdg)== 211) hRcPhiPi->Fill(phi);
      if (abs(pdg)==2212) hRcPhiPr->Fill(phi);
      if (abs(pdg)== 211) hNumberOfMatchedPi->Fill(pt,vNumberOfMatched[ev][ip]);
      if (abs(pdg)== 211) hNumberOfGoodPi->Fill(pt,1);
    }
  }

  new TCanvas;
  auto hEffPtPi = (TH1D*) hRcPtPi->Clone("hEffPtPi");
  hEffPtPi->Divide(hRcPtPi, hSeedPtPi, 1, 1, "B");
  hEffPtPi->Draw();
  double nBest = hRcPtPi->Integral();
  double nSeeds = hSeedPtPi->Integral();
  printf("nBest/nSeeds=%f\n",nBest/nSeeds);
  new TCanvas;
  int nRc = hNumberOfMatchedPi->Integral();
  printf("nRc/nBest=%f\n",nRc/nBest);
  hNumberOfMatchedPi->Divide(hNumberOfMatchedPi,hNumberOfGoodPi,1,1,"B");
  hNumberOfMatchedPi->Draw();
  new TCanvas;
  hNLayers->Draw();
}
