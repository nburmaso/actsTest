R__ADD_INCLUDE_PATH(/home/ekryshen/mpd/actsTest/build/stage/include)
R__ADD_LIBRARY_PATH(/home/ekryshen/mpd/actsTest/build/stage/lib)
R__LOAD_LIBRARY(libactsTestLib.so)
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "ActsFatras/EventData/Barcode.hpp"
#include "tree_summary.C"
#include "MyFtdGeo.h"
#include "MyFtdDetector.h"


const int nStations = 5;
const int nLayersPerStation = 10;
const int shift = 2;
const int minMeasPerCand = 3;
const float minMajFrac = 0.74;
MyFtdGeo* ftdGeo = nullptr;

bool isGoodFtd(int64_t layerMask, int minHits = 5){
  vector<int> nHits(nStations,0); // number of hits per station
  vector<int> nType5(nStations,0); // number of type5 hits per station
  vector<int> nType6(nStations,0); // number of type6 hits per station
  for (int st=0;st<nStations;st++){
    for (int l=0;l<nLayersPerStation;l++) {
      int layerIndex = nLayersPerStation*st+l;
      int type = ftdGeo->GetLayerType(layerIndex);
      if (type==2) continue;
      nHits[st] += TESTBIT(layerMask,layerIndex);
      if (type==5) nType5[st] += TESTBIT(layerMask,layerIndex);
      if (type==6) nType6[st] += TESTBIT(layerMask,layerIndex);
    }
  }
  // printf("%d %d %d %d %d\n",nHits[0],nHits[1],nHits[2],nHits[3],nHits[4]);
  if (nHits[0]<minMeasPerCand) return 0;
  if (nHits[1]<minMeasPerCand) return 0;
  if (nHits[2]<minMeasPerCand) return 0;
  if (nHits[3]<minMeasPerCand) return 0;
  if (nHits[4]<minMeasPerCand) return 0;
  if (nType5[0]==0) return 0;
  if (nType5[1]==0) return 0;
  if (nType5[2]==0) return 0;
  if (nType5[3]==0) return 0;
  if (nType5[4]==0) return 0;
  if (nType6[0]==0) return 0;
  if (nType6[1]==0) return 0;
  if (nType6[2]==0) return 0;
  if (nType6[3]==0) return 0;
  if (nType6[4]==0) return 0;
  return 1;
}

bool isGoodRecoFtd(int64_t layerMask, int minHits = 5){
  vector<int> nHits(nStations,0); // number of hits per station
  for (int st=0;st<nStations;st++){
    for (int l=0;l<nLayersPerStation;l++) {
      int layerIndex = nLayersPerStation*st+l;
      int type = ftdGeo->GetLayerType(layerIndex);
      if (type==2) continue;
      nHits[st] += TESTBIT(layerMask,layerIndex);    }
  }
  if (nHits[0]+nHits[1]+nHits[2]+nHits[3]+nHits[4]<minHits) return 0;
  if (nHits[0]<1) return 0;
  if (nHits[1]<1) return 0;
  if (nHits[2]<1) return 0;
  if (nHits[3]<1) return 0;
  if (nHits[4]<1) return 0;
  return 1;
}

bool isGoodSeed(int64_t layerMask, int minHits = 5){
  vector<int> nHits(nStations,0); // number of hits per station
  vector<int> nType5(nStations,0); // number of type5 hits per station
  vector<int> nType6(nStations,0); // number of type6 hits per station
  for (int st=0;st<nStations;st++){
    for (int l=0;l<nLayersPerStation;l++) {
      int layerIndex = nLayersPerStation*st+l;
      int type = ftdGeo->GetLayerType(layerIndex);
      if (type==2) continue;
      nHits[st] += TESTBIT(layerMask,layerIndex);
      if (type==5) nType5[st] += TESTBIT(layerMask,layerIndex);
      if (type==6) nType6[st] += TESTBIT(layerMask,layerIndex);
    }
  }
  if (nHits[0]<minMeasPerCand) return 0;
  if (nHits[2]<minMeasPerCand) return 0;
  if (nHits[4]<minMeasPerCand) return 0;
  if (nType5[0]==0) return 0;
  if (nType5[2]==0) return 0;
  if (nType5[4]==0) return 0;
  if (nType6[0]==0) return 0;
  if (nType6[2]==0) return 0;
  if (nType6[4]==0) return 0;
  return 1;
}


void analyse_performance(TString dir = "../build/ruvdup90/", double etaMean = 1.75, double etaDif = 0.2, bool refit = 0, bool trackable = 1){
// gStyle->SetOptStat(0);
  ftdGeo = new MyFtdGeo();
  #define axisPt 20,0.,1.
  #define axisPhi 90,-M_PI,M_PI
  #define axisEta 50,1.5,2.0
  TH2D* hPtResVsPtPi = new TH2D("hPtResVsPtPi","",axisPt,2000,-1.,1.);  
  TH2D* hPtResVsPtPr = new TH2D("hPtResVsPtPr","",axisPt,2000,-1.,1.);  
  TH1D* hSeedPtPi = new TH1D("hSeedPtPi","",axisPt);
  TH1D* hSeedPtPr = new TH1D("hSeedPtPr","",axisPt);
  TH1D* hSeedPhiPi = new TH1D("hSeedPhiPi","",axisPhi);
  TH1D* hSeedPhiPr = new TH1D("hSeedPhiPr","",axisPhi);
  TH1D* hSeedEtaPi = new TH1D("hSeedEtaPi","",axisEta);
  TH1D* hSeedEtaPr = new TH1D("hSeedEtaPr","",axisEta);
  TH1D* hSeedablePtPi = new TH1D("hSeedablePtPi","",axisPt);
  TH1D* hSeedablePtPr = new TH1D("hSeedablePtPr","",axisPt);
  TH1D* hSeedablePhiPi = new TH1D("hSeedablePhiPi","",axisPhi);
  TH1D* hSeedablePhiPr = new TH1D("hSeedablePhiPr","",axisPhi);
  TH1D* hSeedableEtaPi = new TH1D("hSeedableEtaPi","",axisEta);
  TH1D* hSeedableEtaPr = new TH1D("hSeedableEtaPr","",axisEta);
  TH1D* hRcPtPi = new TH1D("hRcPtPi","",axisPt);
  TH1D* hRcPtPr = new TH1D("hRcPtPr","",axisPt);
  TH1D* hRcPhiPi = new TH1D("hRcPhiPi","",axisPhi);
  TH1D* hRcPhiPr = new TH1D("hRcPhiPr","",axisPhi);
  TH1D* hRcEtaPi = new TH1D("hRcEtaPi","",axisEta);
  TH1D* hRcEtaPr = new TH1D("hRcEtaPr","",axisEta);
  TH1D* hNumberOfMatchedPi = new TH1D("hNumberOfMatchedPi","",axisPt);
  TH1D* hLayers = new TH1D("hLayers","",40,0,40);
  TH1D* hNLayers = new TH1D("hNLayers","layers",40,0,40);

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

  vector<vector<int64_t>> vFtdLayerMask(nEvents);
  vector<vector<int>> vMatched(nEvents);
  vector<vector<int>> vNumberOfMatched(nEvents);
  vector<vector<int>> vSeeds(nEvents);
  vector<vector<int>> vSpoints0(nEvents);
  vector<vector<int>> vSpoints2(nEvents);
  vector<vector<int>> vSpoints4(nEvents);
  vector<vector<TVector3>> vRcVecP(nEvents);  
  for (int ev=0;ev<nEvents;ev++){ // unordered events (note: ev is not thread safe)
    tPart->GetEntry(ev);
    int nParts = part_pdg->size();
    // counting particles from 1 => increase allocated vector size by 1
    vFtdLayerMask[part_event_id].resize(nParts+1,0);
    vMatched[part_event_id].resize(nParts+1,0);
    vMatched[part_event_id].resize(nParts+1,0);
    vSpoints0[part_event_id].resize(nParts+1,0);    
    vSpoints2[part_event_id].resize(nParts+1,0);    
    vSpoints4[part_event_id].resize(nParts+1,0);    
    vSeeds[part_event_id].resize(nParts+1,0);
    vNumberOfMatched[part_event_id].resize(nParts+1,0);
    vRcVecP[part_event_id].resize(nParts+1);
  }

  printf("fill match array\n");
  TFile* fTrack = new TFile(TString(dir + (refit ? "trackrefit.root" : "tracksummary.root")) );
  TTree* tTrack = (TTree*) fTrack->Get("tracksummary");
  SetBranchAddresses(tTrack);

  for (int ev=0; ev<nEvents; ev++){ // unordered events (note: ev is not thread safe)
    if (ev%100==0) printf("Event = %d\n",ev);
    tTrack->GetEntry(ev);
    for (int it=0; it<m_majorityParticleId->size(); it++){
      int ip = ActsFatras::Barcode().withData(m_majorityParticleId->at(it)).particle();
      if (ip<1) continue;
      if (vMatched[m_eventNr][ip]<1) vMatched[m_eventNr][ip]=1;
      if (!m_hasFittedParams->at(it)) continue;
      if (vMatched[m_eventNr][ip]<2) vMatched[m_eventNr][ip]=2;
      auto& layers = m_measurementLayer->at(it);
      int64_t trackFtdLayerMask = 0;
      for (int il=0;il<layers.size();il++) SETBIT(trackFtdLayerMask, (layers[il]-shift));
      for (int il=0;il<layers.size();il++) hLayers->Fill(layers[il]);
      if (trackable && !isGoodRecoFtd(trackFtdLayerMask)) continue;
      if (vMatched[m_eventNr][ip]<3) vMatched[m_eventNr][ip]=3;
      vNumberOfMatched[m_eventNr][ip]++;
      double qp = m_eQOP_fit->at(it);
      double theta = m_eTHETA_fit->at(it);
      double phi = m_ePHI_fit->at(it);
      vRcVecP[m_eventNr][ip].SetMagThetaPhi(fabs(1./qp), theta, phi);
      hNLayers->Fill(layers.size());
    }
  }

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

  vector<vector<vector<uint32_t>>> vMeasParticleIds(nEvents);
  printf("fill measurement particle ids + map fired layers per particle\n");
  int previous_event = -1;
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if (meas_event_id!=previous_event) {
      if (meas_event_id%100==0) printf("Event=%d\n", meas_event_id);
      previous_event = meas_event_id;
    }
    vMeasParticleIds[meas_event_id].push_back(vector<uint32_t>(meas_particles.size()));
    for (int i=0;i<meas_particles.size();i++){
      int ip = meas_particles[i][2]; // counted from 1
      vMeasParticleIds[meas_event_id].back()[i] = ip;
      SETBIT(vFtdLayerMask[meas_event_id][ip], (meas_layer_id - shift));
    }
  }

  TFile* fSpacepoints = new TFile(dir + "spacepoints.root");
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
    int layerIndex = Acts::GeometryIdentifier(sgeometry_id).layer() - shift;
    int station = ftdGeo->GetLayerStation(layerIndex);
    int majorityId = trunc(smajority);
    float majFrac = (smajority - majorityId) * 10.;
    printf("majorityId=%d majFrac=%f\n", majorityId, majFrac);
    if (majFrac > minMajFrac) {
      if (station==0) vSpoints0[sevent_id][majorityId] += 1;
      if (station==2) vSpoints2[sevent_id][majorityId] += 1;
      if (station==4) vSpoints4[sevent_id][majorityId] += 1;
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

  printf("loop over seeds\n");
  int previous_seed_event_id = -1;
  for (int is=0;is<tSeeds->GetEntries();is++){
    tSeeds->GetEntry(is);
    if (previous_seed_event_id!=seed_event_id){
      previous_seed_event_id = seed_event_id;
      if (seed_event_id%100==0) printf("%d\n",seed_event_id);
    }
    // TODO use majority id from spacepoints
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

  for (int ev=0;ev<nEvents;ev++){ // unordered events (note: ev is not thread safe)
    tPart->GetEntry(ev);
    for (int i=0;i<part_pdg->size();i++){ // particles
      int ip = i+1;
      if (part_mid->at(i)>0) continue; // only primaries
      int64_t ftdLayerMask = vFtdLayerMask[part_event_id][ip];
      if (trackable && (!isGoodFtd(ftdLayerMask) || !isGoodSeed(ftdLayerMask))) continue;    
      int pdg = part_pdg->at(i);
      float vz = part_vz->at(i);
      float pt = part_pt->at(i);
      float eta = part_eta->at(i);
      float phi = part_phi->at(i);    
      if (abs(eta-etaMean)>etaDif || abs(vz)>1.) continue;
      // TODO write seedable selection based on availability of spacepoints
      if (vSpoints0[part_event_id][ip]==0) continue;
      if (vSpoints2[part_event_id][ip]==0) continue;
      if (vSpoints4[part_event_id][ip]==0) continue;
      printf("here\n");
      if (abs(pdg)== 211) hSeedablePtPi->Fill(pt);
      if (abs(pdg)==2212) hSeedablePtPr->Fill(pt);
      if (abs(pdg)== 211) hSeedablePhiPi->Fill(phi);
      if (abs(pdg)==2212) hSeedablePhiPr->Fill(phi);
      if (abs(pdg)== 211) hSeedableEtaPi->Fill(eta);
      if (abs(pdg)==2212) hSeedableEtaPr->Fill(eta);
      if (vSeeds[part_event_id][ip]==0) continue;
      if (vSeeds[part_event_id][ip]>1) printf("Warning: seeds = %d\n",vSeeds[part_event_id][ip]);
      if (abs(pdg)== 211) hSeedPtPi->Fill(pt);
      if (abs(pdg)==2212) hSeedPtPr->Fill(pt);
      if (abs(pdg)== 211) hSeedPhiPi->Fill(phi);
      if (abs(pdg)==2212) hSeedPhiPr->Fill(phi);
      if (abs(pdg)== 211) hSeedEtaPi->Fill(eta);
      if (abs(pdg)==2212) hSeedEtaPr->Fill(eta);
      if (vMatched[part_event_id][ip]<3) continue;
      if (abs(pdg)== 211) hRcPtPi->Fill(pt);
      if (abs(pdg)==2212) hRcPtPr->Fill(pt);
      if (abs(pdg)== 211) hRcPhiPi->Fill(phi);
      if (abs(pdg)==2212) hRcPhiPr->Fill(phi);
      if (abs(pdg)== 211) hRcEtaPi->Fill(eta);
      if (abs(pdg)==2212) hRcEtaPr->Fill(eta);
      if (abs(pdg)== 211) hNumberOfMatchedPi->Fill(pt,vNumberOfMatched[part_event_id][ip]);
      double ptRC = vRcVecP[part_event_id][ip].Perp();
      if (abs(pdg)== 211) hPtResVsPtPi->Fill(pt,(ptRC-pt)/pt);
      if (abs(pdg)==2212) hPtResVsPtPr->Fill(pt,(ptRC-pt)/pt);
    }
  }

  double nSeedable = hSeedablePtPi->Integral();
  double nSeeds = hSeedPtPi->Integral();
  double nBest = hRcPtPi->Integral();
  double nRc = hNumberOfMatchedPi->Integral();
  printf("nSeeds/nSeedable=%f\n",nSeeds/nSeedable);
  printf("nBest/nSeeds=%f\n",nBest/nSeeds);
  printf("nRc/nBest=%f\n",nRc/nBest);

  new TCanvas;
  hLayers->Draw();

  new TCanvas;
  hNLayers->Draw();

  new TCanvas;
  hNumberOfMatchedPi->Divide(hNumberOfMatchedPi,hRcPtPi,1,1,"B");
  hNumberOfMatchedPi->Draw();

  new TCanvas;
  auto hSeedEffPtPi = (TH1D*) hSeedPtPi->Clone("hSeedEffPtPi");
  hSeedEffPtPi->Divide(hSeedPtPi, hSeedablePtPi, 1, 1, "B");
  hSeedEffPtPi->Draw();

  new TCanvas;
  auto hSeedEffPhiPi = (TH1D*) hSeedPhiPi->Clone("hSeedEffPhiPi");
  hSeedEffPhiPi->Divide(hSeedPhiPi, hSeedablePhiPi, 1, 1, "B");
  hSeedEffPhiPi->SetMinimum(0);
  hSeedEffPhiPi->Draw();

  new TCanvas;
  auto hSeedEffEtaPi = (TH1D*) hSeedEtaPi->Clone("hSeedEffEtaPi");
  hSeedEffEtaPi->Divide(hSeedEtaPi, hSeedableEtaPi, 1, 1, "B");
  hSeedEffEtaPi->Draw();

  new TCanvas;
  auto hRcEffPtPi = (TH1D*) hRcPtPi->Clone("hRcEffPtPi");
  hRcEffPtPi->Divide(hRcPtPi, hSeedPtPi, 1, 1, "B");
  hRcEffPtPi->Draw();

  new TCanvas;
  auto hRcEffPhiPi = (TH1D*) hRcPhiPi->Clone("hRcEffPhiPi");
  hRcEffPhiPi->Divide(hRcPhiPi, hSeedPhiPi, 1, 1, "B");
  hRcEffPhiPi->SetMinimum(0);
  hRcEffPhiPi->Draw();

  new TCanvas;
  auto hRcEffEtaPi = (TH1D*) hRcEtaPi->Clone("hRcEffEtaPi");
  hRcEffEtaPi->Divide(hRcEtaPi, hSeedEtaPi, 1, 1, "B");
  hRcEffEtaPi->Draw();

  auto* fout = new TFile(dir +Form("tracking_performance_%.2f.root",etaMean), "recreate");
  hSeedPtPi->Write();
  hSeedPtPr->Write();
  hSeedPhiPi->Write();
  hSeedPhiPr->Write();
  hSeedEtaPi->Write();
  hSeedEtaPr->Write();
  hSeedablePtPi->Write();
  hSeedablePtPr->Write();
  hSeedablePhiPi->Write();
  hSeedablePhiPr->Write();
  hSeedableEtaPi->Write();
  hSeedableEtaPr->Write();
  hRcPtPi->Write();
  hRcPtPr->Write();
  hRcPhiPi->Write();
  hRcPhiPr->Write();
  hRcEtaPi->Write();
  hRcEtaPr->Write();
  hNumberOfMatchedPi->Write();
  hLayers->Write();
  hNLayers->Write();
  hPtResVsPtPi->Write();
  hPtResVsPtPr->Write();

  fout->Close();

}
