#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "vector"
#include "style.h"
#include "tree_summary.C"
#include "tree_trackstates.C"


bool isGoodFtd(int64_t layerMask, int minHits = 5, int shift = 3){
  int nHits[5] = {0,0,0,0,0}; // number of hits per station
  for (int st=0;st<5;st++){
    nHits[st] += (layerMask & (1ull << (7*st + shift + 0))) > 0;
    nHits[st] += (layerMask & (1ull << (7*st + shift + 1))) > 0;
    nHits[st] += (layerMask & (1ull << (7*st + shift + 2))) > 0;
    nHits[st] += (layerMask & (1ull << (7*st + shift + 4))) > 0;
    nHits[st] += (layerMask & (1ull << (7*st + shift + 5))) > 0;
    nHits[st] += (layerMask & (1ull << (7*st + shift + 6))) > 0;
  }
  // printf("%d %d %d %d %d\n",nHits[0],nHits[1],nHits[2],nHits[3],nHits[4]);
  if (nHits[0]<2) return 0;
  if (nHits[1]<2) return 0;
  if (nHits[2]<2) return 0;
  if (nHits[3]<2) return 0;
  if (nHits[4]<2) return 0;
  return 1;
}

bool isGoodSeed(int64_t layerMask, int minHits = 5, int shift = 3){
  int nSeeds[5] = {0,0,0,0,0}; // number of seeds per station
  int nSeedsB = (layerMask & (1ull << (shift -1))) > 0;
  int nSeedsF = (layerMask & (1ull << (7*4 + shift + 7))) > 0;
  for (int st=0;st<5;st++){
    nSeeds[st] += (layerMask & (1ull << (7*st + shift + 3))) > 0;
  }
  if (nSeedsB<1) return 0;
  if (nSeedsF<1) return 0;
  // if (nSeeds[0]<1) return 0;
  if (nSeeds[2]<1) return 0;
  // if (nSeeds[4]<1) return 0;
  return 1;
}

bool isGoodRecoFtd(int64_t layerMask, int minHits = 5, int shift = 3){
  int nHits[5] = {0,0,0,0,0}; // number of hits per station
  for (int st=0;st<5;st++){
    nHits[st] += (layerMask & (1ull << (7*st + shift + 0))) > 0;
    nHits[st] += (layerMask & (1ull << (7*st + shift + 1))) > 0;
    nHits[st] += (layerMask & (1ull << (7*st + shift + 2))) > 0;
    nHits[st] += (layerMask & (1ull << (7*st + shift + 4))) > 0;
    nHits[st] += (layerMask & (1ull << (7*st + shift + 5))) > 0;
    nHits[st] += (layerMask & (1ull << (7*st + shift + 6))) > 0;
  }
  if (nHits[0]<1) return 0;
  if (nHits[1]<1) return 0;
  if (nHits[2]<1) return 0;
  if (nHits[3]<1) return 0;
  if (nHits[4]<1) return 0;
  return 1;
}

void analyse_tracks(TString dir = "../build/test/", double etaMean = 1.7, double etaDif = 0.1, bool refit = 0, bool trackable = 1){
  TFile* fPart = new TFile(dir + "particles.root");
  TTree* tPart = (TTree*) fPart->Get("particles");
  uint32_t part_event_id;
  auto part_id  = new std::vector<unsigned long>; 
  auto part_pdg = new vector<int>;
  auto part_vz  = new std::vector<float>;
  auto part_pt  = new std::vector<float>;
  auto part_eta = new std::vector<float>;
  auto part_phi = new std::vector<float>;
  tPart->SetBranchAddress("event_id",&part_event_id);
  tPart->SetBranchAddress("particle",&part_id);  
  tPart->SetBranchAddress("particle_type",&part_pdg);
  tPart->SetBranchAddress("vz",&part_vz);
  tPart->SetBranchAddress("pt",&part_pt);  
  tPart->SetBranchAddress("eta",&part_eta);  
  tPart->SetBranchAddress("phi",&part_phi);  
  int nEvents = tPart->GetEntries();

  TFile* fTracks = new TFile(dir + "tracks.root");
  TTree* tTracks = (TTree*) fTracks->Get("tracks");
  std::uint32_t track_eventNr{0};
  std::vector<std::uint32_t> track_trackNr;       auto ptrack_trackNr = &track_trackNr;
  std::vector<float> track_chi2Sum;               auto ptrack_chi2Sum = &track_chi2Sum;
  std::vector<unsigned int> track_NDF;            auto ptrack_NDF = &track_NDF;
  std::vector<float> track_eLOC0;                 auto ptrack_eLOC0 = &track_eLOC0;
  std::vector<float> track_eLOC1;                 auto ptrack_eLOC1 = &track_eLOC1;
  std::vector<float> track_ePHI;                  auto ptrack_ePHI = &track_ePHI;
  std::vector<float> track_eTHETA;                auto ptrack_eTHETA = &track_eTHETA;
  std::vector<float> track_eQOP;                  auto ptrack_eQOP = &track_eQOP;
  std::vector<std::vector<int>> track_measurementIds; auto ptrack_measurementIds = &track_measurementIds;
  
  tTracks->SetBranchAddress("event_nr", &track_eventNr);
  tTracks->SetBranchAddress("track_nr", &ptrack_trackNr);
  tTracks->SetBranchAddress("chi2Sum", &ptrack_chi2Sum);
  tTracks->SetBranchAddress("NDF", &ptrack_NDF);
  tTracks->SetBranchAddress("eLOC0_fit", &ptrack_eLOC0);
  tTracks->SetBranchAddress("eLOC1_fit", &ptrack_eLOC1);
  tTracks->SetBranchAddress("ePHI_fit", &ptrack_ePHI);
  tTracks->SetBranchAddress("eTHETA_fit", &ptrack_eTHETA);
  tTracks->SetBranchAddress("eQOP_fit", &ptrack_eQOP);
  tTracks->SetBranchAddress("measurementIds", &ptrack_measurementIds);

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

  printf("fill measurement particle ids + map fired layers per particle\n");

  vector<map<int,int64_t>> vEventParticleFtdLayerMask(nEvents);
  map<int,int> mapFirstMeaurementIndex;
  vector<map<int,int>> mapMeasParticles(nEvents);
  vector<map<int,int>> mapMeasLayers(nEvents);
  int previous_event_nr = -1;
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if (meas_event_id!=previous_event_nr) {
      if (meas_event_id%10==0) printf("Event=%d\n", meas_event_id);
      previous_event_nr = meas_event_id;
      mapFirstMeaurementIndex[meas_event_id] = im;
    }
    int ip = meas_particles[0][2] - 1;
    int imRel = im - mapFirstMeaurementIndex[meas_event_id];
    mapMeasParticles[meas_event_id][imRel] = ip;
    mapMeasLayers[meas_event_id][imRel] = meas_layer_id;
    auto& mapParticleFtdLayerMask = vEventParticleFtdLayerMask[meas_event_id];
    if (auto m = mapParticleFtdLayerMask.find(ip); m == mapParticleFtdLayerMask.end()) mapParticleFtdLayerMask[ip] = 0;
    mapParticleFtdLayerMask[ip] |= (1ull << meas_layer_id);
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

  vector<vector<int>> vSeeds(nEvents);
  vector<vector<int>> vMatched(nEvents);
  vector<vector<int>> vNumberOfMatched(nEvents);

  for (int ev=0; ev<nEvents; ev++){
    tPart->GetEntry(ev);
    auto nParticles = part_id->size();
    vSeeds[ev].assign(nParticles,0);
    vMatched[ev].assign(nParticles,0);
    vNumberOfMatched[ev].assign(nParticles,0);
  }

  printf("loop over seeds\n");
  int previous_seed_event_id = -1;
  for (int is=0;is<tSeeds->GetEntries();is++){
    tSeeds->GetEntry(is);
    if (previous_seed_event_id!=seed_event_id){
      previous_seed_event_id = seed_event_id;
      if (seed_event_id%1==0) printf("Seed event %d\n",seed_event_id);
    }
    auto& partIds = mapMeasParticles[seed_event_id];
    auto id1 = partIds[measId1];
    auto id2 = partIds[measId2];
    auto id3 = partIds[measId3];
    //printf("%d %d %d %d\n",seed_event_id, id1, id2, id3);
    if (id1==id2 && id2==id3) vSeeds[seed_event_id][id1]+=1;
  }

  TH1D* hRcPtPi = new TH1D("hRcPtPi","",100,0.,1.);
  TH1D* hRcPtPr = new TH1D("hRcPtPr","",100,0.,1.);
  TH1D* hRcPhiPi = new TH1D("hRcPhiPi","",180,-M_PI,M_PI);
  TH1D* hRcPhiPr = new TH1D("hRcPhiPr","",180,-M_PI,M_PI);
  TH1D* hNumberOfMatched = new TH1D("hNumberOfMatched","",100,0.,1.);
  TH1D* hNumberOfGood = new TH1D("hNumberOfGood","",100,0.,1.);

  printf("loop over tracks\n");
  for (int ev=0; ev<tTracks->GetEntries(); ev++){
    if (ev%1==0) printf("Event=%d\n", ev);
    tTracks->GetEntry(ev);
    tPart->GetEntry(ev);

    int nTracks = track_trackNr.size();
    map<int,set<int>> tracksPerMeasurement;
    vector<float> trackChi2(nTracks,0);
    vector<int64_t> trackFtdLayerMask(nTracks,0);
    vector<int> majorityParticleId(nTracks,0);
    vector<int> majorityHits(nTracks,0);
    for (int it=0;it<nTracks;it++){
      //printf("%d %f\n", track_trackNr[it], track_chi2Sum[it]/track_NDF[it]);
      trackChi2[it] = track_NDF[it]>0 ? track_chi2Sum[it]/track_NDF[it] : 1e10;
      map<int,int> mapParticles;
      for (auto measId : track_measurementIds[it]) {
        printf("  %d %d\n", measId, mapMeasParticles[ev][measId]);
        tracksPerMeasurement[measId].insert(it);
        trackFtdLayerMask[it] |= (1ull << mapMeasLayers[ev][measId]);
        int partId = mapMeasParticles[ev][measId];
        mapParticles[partId]++;
      }
      // find majority partId
      int majorityPartId = -1;
      int nMaxHits = 0;
      for (auto part : mapParticles){
        if (part.second <= nMaxHits) continue;
        majorityPartId = part.first;
        nMaxHits = part.second;
      }
      majorityParticleId[it] = majorityPartId;
      majorityHits[it] = nMaxHits;
    }
    
    // printf("tracksPerMeasurement: \n");
    // for (auto& meas: tracksPerMeasurement) {
    //   printf("%d %d\n",meas.first, meas.second.size());
    // }

    vector<int> sharedMeasurementsPerTrack(nTracks, 0);
    set<int> selectedTracks;
    for (int it=0;it<nTracks;it++){
      for (auto measId : track_measurementIds[it]) {
        if (tracksPerMeasurement[measId].size() > 1) {
          sharedMeasurementsPerTrack[it]++;
        }
      }
      printf("Track %d, meas: %d, shared: %d, chi2=%f\n",it, track_measurementIds[it].size(), sharedMeasurementsPerTrack[it], trackChi2[it]);
      selectedTracks.insert(it);
    }

    auto sharedMeasurementsComperator = [&sharedMeasurementsPerTrack](std::size_t a, std::size_t b) {
      return sharedMeasurementsPerTrack[a] < sharedMeasurementsPerTrack[b];
    };

    auto trackComperator = [&sharedMeasurementsPerTrack,&track_measurementIds,&trackChi2](std::size_t a, std::size_t b) {
      auto relSharedA = 1.*sharedMeasurementsPerTrack[a]/track_measurementIds[a].size();
      auto relSharedB = 1.*sharedMeasurementsPerTrack[b]/track_measurementIds[b].size();
      if (relSharedA != relSharedB) return relSharedA < relSharedB;
      if (track_measurementIds[a].size() == track_measurementIds[b].size()) return trackChi2[a] < trackChi2[b];
      return track_measurementIds[a].size() > track_measurementIds[b].size();
    };

    int maximumSharedHits = 1;
    for (int i = 0; i < 1000; i++) {
      if (selectedTracks.empty()) break;
      auto maximumSharedMeasurements = *std::max_element(selectedTracks.begin(), selectedTracks.end(), sharedMeasurementsComperator);
      printf("track %d with maximumSharedMeasurements=%d\n",maximumSharedMeasurements,sharedMeasurementsPerTrack[maximumSharedMeasurements]);
      if (sharedMeasurementsPerTrack[maximumSharedMeasurements] <= maximumSharedHits) break;
      auto badTrack = *std::max_element(selectedTracks.begin(), selectedTracks.end(), trackComperator);
      printf("badTrack=%d %zu\n",badTrack, track_measurementIds[badTrack].size());
      // remove bad track
      for (auto im : track_measurementIds[badTrack]) {
        tracksPerMeasurement[im].erase(badTrack);
        if (tracksPerMeasurement[im].size()>1) continue;
        auto it = *tracksPerMeasurement[im].begin();
        sharedMeasurementsPerTrack[it]--;
      }
      selectedTracks.erase(badTrack);
    }
    
    vector<int> isSelected(nTracks,0);
    for (auto it : selectedTracks) isSelected[it] = 1;

    int nParticles = part_id->size();

    for (int it=0;it<nTracks;it++){
      if (!isSelected[it]) continue;
      if (track_NDF[it]==0) continue;
      printf("track=%d chi2=%f majPartId=%d majHits=%d\n", track_trackNr[it], track_chi2Sum[it]/track_NDF[it], majorityParticleId[it],majorityHits[it]);
      auto ip = majorityParticleId[it];
      if (vMatched[ev][ip]<2) vMatched[ev][ip]=2;
      if (!isGoodRecoFtd(trackFtdLayerMask[it])) continue;
      if (vMatched[ev][ip]<3) vMatched[ev][ip]=3;
      vNumberOfMatched[ev][ip]++;
    }

    for (int ip=0;ip<nParticles;ip++){ // particles
      // if (part_mid->at(ip)>0) continue; // only primaries
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
      hNumberOfMatched->Fill(pt,vNumberOfMatched[ev][ip]);
      hNumberOfGood->Fill(pt,1);
    }
  } // events

  new TCanvas;
  auto hEffPtPi = (TH1D*) hRcPtPi->Clone("hEffPtPi");
  hEffPtPi->Divide(hRcPtPi, hSeedPtPi, 1, 1, "B");
  hEffPtPi->Draw();
  double nBest = hRcPtPi->Integral();
  double nSeeds = hSeedPtPi->Integral();
  printf("nBest/nSeeds=%f\n",nBest/nSeeds);
  new TCanvas;
  int nRc = hNumberOfMatched->Integral();
  printf("nRc/nBest=%f\n",nRc/nBest);
  hNumberOfMatched->Divide(hNumberOfMatched,hNumberOfGood,1,1,"B");
  hNumberOfMatched->Draw();

}