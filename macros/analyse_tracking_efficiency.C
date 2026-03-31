#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "ActsFatras/EventData/Barcode.hpp"
#include "tree_summary.C"


const int nStations = 5;
const int nLayersPerStation = 9;


bool isGoodFtd(int64_t layerMask, int minHits = 5){
  vector<int> nHits(nStations,0); // number of hits per station
  for (int st=0;st<nStations;st++){
    for (int l=0;l<nLayersPerStation;l++) {
      if (l==nLayersPerStation/2) continue;
      nHits[st] += ((layerMask & (1ull << (nLayersPerStation*st+l))) > 0);
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
  vector<int> nHits(nStations,0); // number of hits per station
  for (int st=0;st<nStations;st++){
    for (int l=0;l<nLayersPerStation;l++) {
      if (l==nLayersPerStation/2) continue;
      nHits[st] += ((layerMask & (1ull << (nLayersPerStation*st+l))) > 0);
    }
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
  vector<int> nSeeds(nStations,0); // number of hits per station
  for (int st=0;st<nStations;st++){
    for (int l=0;l<nLayersPerStation;l++) {
      if (l==nLayersPerStation/2) continue;
      nSeeds[st] += ((layerMask & (1ull << (nLayersPerStation*st+l))) > 0);
    }
  }
  if (nSeeds[0]<3) return 0;
  if (nSeeds[2]<3) return 0;
  if (nSeeds[4]<3) return 0;
  return 1;
}

void analyse_tracking_efficiency(TString dir = "../build/test", double etaMean = 1.9, double etaDif = 0.05, bool refit = 0, bool trackable = 1, int nPrimaries = 1){
  dir.Append("/");
  // int shift = 3; //  isroc = 1;  isframe = 1; with fake pre-layer
  int shift = 2; //  isroc = 1;  isframe = 1; without fake pre-layer

  // setup particles
  TFile* fPart = new TFile(TString(dir + "particles.root"));
  TTree* tPart = (TTree*) fPart->Get("particles");
  UInt_t part_event_id;
  auto part_id  = new std::vector<uint32_t>; 
  auto part_pdg = new vector<int>;
  auto part_vz  = new std::vector<float>;
  auto part_px  = new std::vector<float>;
  auto part_py  = new std::vector<float>;
  auto part_pz  = new std::vector<float>;
  auto part_pt  = new std::vector<float>;
  auto part_eta = new std::vector<float>;
  auto part_phi = new std::vector<float>;
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

  int nEvents = tPart->GetEntries();

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
    if (ip < 0 || ip>=nPrimaries) continue;
    int64_t ftdLayerMask = vEventParticleFtdLayerMask[seed_event_id][ip];
    if (trackable && (!isGoodFtd(ftdLayerMask) || !isGoodSeed(ftdLayerMask))) continue;    
    int pdg = part_pdg->at(ip);
    float vz = part_vz->at(ip);
    float pt = part_pt->at(ip);
    float eta = part_eta->at(ip);
    float phi = part_phi->at(ip);    
    if (abs(eta-etaMean)<etaDif && abs(vz)<1.) {
      if (abs(pdg)== 211) hSeedPtPi->Fill(pt);
      if (abs(pdg)==2212) hSeedPtPr->Fill(pt);
      if (abs(pdg)== 211) hSeedPhiPi->Fill(phi);
      if (abs(pdg)==2212) hSeedPhiPr->Fill(phi);
    }
  }


  TFile* fTrack = new TFile(TString(dir + (refit ? "trackrefit.root" : "tracksummary.root")) );
  TTree* tTrack = (TTree*) fTrack->Get("tracksummary");
  SetBranchAddresses(tTrack);

  TH1D* hChi2all = new TH1D("hChi2all","",1000,0.,100.);
  TH1D* hChi2sel = new TH1D("hChi2sel","",1000,0.,100.);

  TH1D* hMcPtPi = new TH1D("hMcPtPi","",100,0.,1.);
  TH1D* hMcPtPr = new TH1D("hMcPtPr","",100,0.,1.);
  TH1D* hRcPtPi = new TH1D("hRcPtPi","",100,0.,1.);
  TH1D* hRcPtPr = new TH1D("hRcPtPr","",100,0.,1.);

  TH1D* hMcPhiPi = new TH1D("hMcPhiPi","",180,-M_PI,M_PI);
  TH1D* hMcPhiPr = new TH1D("hMcPhiPr","",180,-M_PI,M_PI);
  TH1D* hRcPhiPi = new TH1D("hRcPhiPi","",180,-M_PI,M_PI);
  TH1D* hRcPhiPr = new TH1D("hRcPhiPr","",180,-M_PI,M_PI);

  TH2D* hPResVsPtPi = new TH2D("hPResVsPtPi","",20,0.,1.,2000,-1.,1.);
  TH2D* hPResVsPtPr = new TH2D("hPResVsPtPr","",20,0.,1.,2000,-1.,1.);

  TH2D* hPtResVsPtPi = new TH2D("hPtResVsPtPi","",20,0.,1.,2000,-1.,1.);
  TH2D* hPtResVsPtPr = new TH2D("hPtResVsPtPr","",20,0.,1.,2000,-1.,1.);

  TH2D* hResQOPvsPtPi  = new TH2D("hResQOPvsPtPi","",20,0.,1.,200,-1,1);
  TH2D* hResDvsPtPi    = new TH2D("hResDvsPtPi","",20,0.,1.,200,-9,9);
  TH2D* hResZvsPtPi    = new TH2D("hResZvsPtPi","",20,0.,1.,200,-9,9);
  TH2D* hPullQOPvsPtPi = new TH2D("hPullQOPvsPtPi","",20,0.,1.,200,-10,10);
  TH2D* hPullDvsPtPi   = new TH2D("hPullDvsPtPi","",20,0.,1.,200,-10,10);
  TH2D* hPullZvsPtPi   = new TH2D("hPullZvsPtPi","",20,0.,1.,200,-10,10);

  TH2D* hResQOPvsPtPr  = new TH2D("hResQOPvsPtPr","",20,0.,1.,200,-1,1);
  TH2D* hResDvsPtPr    = new TH2D("hResDvsPtPr","",20,0.,1.,200,-9,9);
  TH2D* hResZvsPtPr    = new TH2D("hResZvsPtPr","",20,0.,1.,200,-9,9);
  TH2D* hPullQOPvsPtPr = new TH2D("hPullQOPvsPtPr","",20,0.,1.,200,-10,10);
  TH2D* hPullDvsPtPr   = new TH2D("hPullDvsPtPr","",20,0.,1.,200,-10,10);
  TH2D* hPullZvsPtPr   = new TH2D("hPullZvsPtPr","",20,0.,1.,200,-10,10);
  TH1D* hMeasurements  = new TH1D("hMeasurements","",40,0,40.);
  TH1D* hMajorityHits  = new TH1D("hMajorityHits","",40,0,40.);

  TVector3 v;
  TVector3 vMC;
  printf("loop over tracks\n");
  for (int ev=0; ev<tTrack->GetEntries(); ev++){
    if (ev%100==0) printf("Event = %d\n",ev);
    tTrack->GetEntry(ev);
    tPart->GetEntry(ev);
    auto& mapParticleFtdLayerMask = vEventParticleFtdLayerMask[ev];
    for (int ip = 0; ip<part_pdg->size(); ip++){
      int64_t ftdLayerMask = mapParticleFtdLayerMask[ip];
      if (trackable && (!isGoodFtd(ftdLayerMask) || !isGoodSeed(ftdLayerMask))) continue;    
      if (ip < 0 || ip>=nPrimaries) continue;
      int pdg = part_pdg->at(ip);
      float vz = part_vz->at(ip);
      float pt = part_pt->at(ip);
      float eta = part_eta->at(ip);
      float phi = part_phi->at(ip);    
      if (abs(eta-etaMean)<etaDif && abs(vz)<1.) {
        if (abs(pdg)== 211) hMcPtPi->Fill(pt);
        if (abs(pdg)==2212) hMcPtPr->Fill(pt);
        if (abs(pdg)== 211) hMcPhiPi->Fill(phi);
        if (abs(pdg)==2212) hMcPhiPr->Fill(phi);
      }
    } // particle loop

    for (int it=0; it<m_majorityParticleId->size(); it++){
      if (!m_hasFittedParams->at(it)) continue;
      // check measurements and true hits here
      


      hChi2all->Fill(m_chi2Sum->at(it)/m_nMeasurements->at(it));
      auto& layers = m_measurementLayer->at(it);
      int64_t trackFtdLayerMask = 0;
      for (int il=0;il<layers.size();il++) trackFtdLayerMask |= (1ull << (layers[il] - shift));
      if (trackable && !isGoodFtd(trackFtdLayerMask)) continue;
      hMeasurements->Fill(m_nMeasurements->at(it));
      hMajorityHits->Fill(m_nMajorityHits->at(it));
      hChi2sel->Fill(m_chi2Sum->at(it)/m_nMeasurements->at(it));
      double qp = m_eQOP_fit->at(it);
      double theta = m_eTHETA_fit->at(it);
      double phi = m_ePHI_fit->at(it);
      v.SetMagThetaPhi(1./qp, theta, phi);
      double pRC = v.Mag();
      double ptRC = v.Pt();
      auto barcode = ActsFatras::Barcode().withData(m_majorityParticleId->at(it));
      int ip = barcode.particle()-1;
      if (ip < 0 || ip>=nPrimaries) continue;
      int64_t ftdLayerMask = mapParticleFtdLayerMask[ip];
      if (trackable && (!isGoodFtd(ftdLayerMask) || !isGoodSeed(ftdLayerMask))) continue;    
      int pdg = part_pdg->at(ip);
      float vz = part_vz->at(ip);
      float ptMC = part_pt->at(ip);
      float etaMC = part_eta->at(ip);
      float phiMC = part_phi->at(ip);
      vMC.SetPtEtaPhi(ptMC, etaMC, phiMC);
      double pMC = vMC.Mag();
      float qpMC = 1/vMC.Mag();
      if (abs(etaMC-etaMean)<etaDif && abs(vz)<1.) {
        if (abs(pdg)== 211) {
          hRcPtPi->Fill(ptMC);
          hRcPhiPi->Fill(phiMC);
          hPResVsPtPi->Fill(ptMC,(pRC-pMC)/pMC);
          hPtResVsPtPi->Fill(ptMC,(ptRC-ptMC)/ptMC);
          hResQOPvsPtPi->Fill(ptMC,m_eQOP_fit->at(it)-qpMC);
          hResDvsPtPi->Fill(ptMC,m_eLOC0_fit->at(it)/10.);
          hResZvsPtPi->Fill(ptMC,m_eLOC1_fit->at(it)/10.);
          hPullQOPvsPtPi->Fill(ptMC,m_pull_eQOP_fit->at(it));
          hPullDvsPtPi->Fill(ptMC,m_pull_eLOC0_fit->at(it));
          hPullZvsPtPi->Fill(ptMC,m_pull_eLOC1_fit->at(it));
        }
        if (abs(pdg)==2212) {
          hRcPtPr->Fill(ptMC);
          hRcPhiPr->Fill(phiMC);
          hPResVsPtPr->Fill(ptMC,(pRC-pMC)/pMC);
          hPtResVsPtPr->Fill(ptMC,(ptRC-ptMC)/ptMC);
          hResQOPvsPtPr->Fill(ptMC,m_eQOP_fit->at(it)-qpMC);
          hResDvsPtPr->Fill(ptMC,m_eLOC0_fit->at(it)/10.);
          hResZvsPtPr->Fill(ptMC,m_eLOC1_fit->at(it)/10.);
          hPullQOPvsPtPr->Fill(ptMC,m_pull_eQOP_fit->at(it));
          hPullDvsPtPr->Fill(ptMC,m_pull_eLOC0_fit->at(it));
          hPullZvsPtPr->Fill(ptMC,m_pull_eLOC1_fit->at(it));
        }
      } // track loop
    } // event loop
  }
  new TCanvas;
  hMeasurements->Draw();

  new TCanvas;
  hMcPtPi->Draw();

  new TCanvas;
  auto hEffPtPi = (TH1D*) hRcPtPi->Clone("hEffPtPi");
  hEffPtPi->Divide(hRcPtPi, hMcPtPi, 1, 1, "B");
  hEffPtPi->Draw();
  
  new TCanvas;
  hPtResVsPtPi->Draw();

  if (0) {
  new TCanvas;
  hMcPtPr->Draw();

  new TCanvas;
  auto hEffPtPr = (TH1D*) hRcPtPr->Clone("hEffPtPr");
  hEffPtPr->Divide(hRcPtPr, hMcPtPr, 1, 1, "B");
  hEffPtPr->Draw();

  new TCanvas;
  hPtResVsPtPr->Draw();
  }

  new TCanvas;
  hChi2all->SetLineColor(kBlue);
  hChi2sel->SetLineColor(kRed);
  hChi2all->Draw();
  hChi2sel->Draw("same");

  TFile* f = new TFile(TString(dir + (refit ? "tracking_efficiency_refit.root" : "tracking_efficiency.root") ),"update");
  hMajorityHits->Write(Form("hMajorityHits%.0f",etaMean*10));
  hMeasurements->Write(Form("hMeasurements%.0f",etaMean*10));
  hMcPtPi->Write(Form("hMcPtPi%.0f",etaMean*10));
  hMcPtPr->Write(Form("hMcPtPr%.0f",etaMean*10));
  hRcPtPi->Write(Form("hRcPtPi%.0f",etaMean*10));
  hRcPtPr->Write(Form("hRcPtPr%.0f",etaMean*10));
  hMcPhiPi->Write(Form("hMcPhiPi%.0f",etaMean*10));
  hMcPhiPr->Write(Form("hMcPhiPr%.0f",etaMean*10));
  hRcPhiPi->Write(Form("hRcPhiPi%.0f",etaMean*10));
  hRcPhiPr->Write(Form("hRcPhiPr%.0f",etaMean*10));
  hSeedPtPi->Write(Form("hSeedPtPi%.0f",etaMean*10));
  hSeedPtPr->Write(Form("hSeedPtPr%.0f",etaMean*10));
  hSeedPhiPi->Write(Form("hSeedPhiPi%.0f",etaMean*10));
  hSeedPhiPr->Write(Form("hSeedPhiPr%.0f",etaMean*10));
  hPResVsPtPi->Write(Form("hPResVsPtPi%.0f",etaMean*10));
  hPResVsPtPr->Write(Form("hPResVsPtPr%.0f",etaMean*10));
  hPtResVsPtPi->Write(Form("hPtResVsPtPi%.0f",etaMean*10));
  hPtResVsPtPr->Write(Form("hPtResVsPtPr%.0f",etaMean*10));
  hResQOPvsPtPi->Write(Form("hResQOPvsPtPi%.0f",etaMean*10));
  hResDvsPtPi->Write(Form("hResDvsPtPi%.0f",etaMean*10));
  hResZvsPtPi->Write(Form("hResZvsPtPi%.0f",etaMean*10));
  hPullQOPvsPtPi->Write(Form("hPullQOPvsPtPi%.0f",etaMean*10));
  hPullDvsPtPi->Write(Form("hPullDvsPtPi%.0f",etaMean*10));
  hPullZvsPtPi->Write(Form("hPullZvsPtPi%.0f",etaMean*10));
  hResQOPvsPtPr->Write(Form("hResQOPvsPtPr%.0f",etaMean*10));
  hResDvsPtPr->Write(Form("hResDvsPtPr%.0f",etaMean*10));
  hResZvsPtPr->Write(Form("hResZvsPtPr%.0f",etaMean*10));
  hPullQOPvsPtPr->Write(Form("hPullQOPvsPtPr%.0f",etaMean*10));
  hPullDvsPtPr->Write(Form("hPullDvsPtPr%.0f",etaMean*10));
  hPullZvsPtPr->Write(Form("hPullZvsPtPr%.0f",etaMean*10));
  f->Close();
}
