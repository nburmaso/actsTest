#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "ActsFatras/EventData/Barcode.hpp"


void analyse_seed_efficiency(TString dir = "mpd", double etaMean = 1.6, double etaDif = 0.1){
  dir.Append("/");
  // setup particles
  TFile* fPart = new TFile(TString(dir + "particles.root"));
  TTree* tPart = (TTree*) fPart->Get("particles");
  UInt_t part_event_id;
  auto part_id  = new std::vector<unsigned long>; 
  auto part_pdg = new vector<int>;
  auto part_vz  = new std::vector<float>;
  auto part_px  = new std::vector<float>;
  auto part_py  = new std::vector<float>;
  auto part_pz  = new std::vector<float>;
  auto part_pt  = new std::vector<float>;
  auto part_eta = new std::vector<float>;
  tPart->SetBranchAddress("event_id",&part_event_id);
  tPart->SetBranchAddress("particle_id",&part_id);
  tPart->SetBranchAddress("particle_type",&part_pdg);
  tPart->SetBranchAddress("vz",&part_vz);
  tPart->SetBranchAddress("px",&part_px);
  tPart->SetBranchAddress("py",&part_py);
  tPart->SetBranchAddress("pz",&part_pz);  
  tPart->SetBranchAddress("pt",&part_pt);  
  tPart->SetBranchAddress("eta",&part_eta);  

  // setup seeds
  TFile* fSeeds = new TFile(TString(dir + "seeds.root"));
  TTree* tSeeds = (TTree*) fSeeds->Get("seeds");
  UInt_t seed_event_id;
  ULong64_t meas_id_1;
  ULong64_t meas_id_2;
  ULong64_t meas_id_3;
  tSeeds->SetBranchAddress("event_id", &seed_event_id);
  tSeeds->SetBranchAddress("measurement_id_1", &meas_id_1);
  tSeeds->SetBranchAddress("measurement_id_2", &meas_id_2);
  tSeeds->SetBranchAddress("measurement_id_3", &meas_id_3);

  // setup measurements
  TFile* fMeas = new TFile(TString(dir + "measurements.root"));
  TTree* tMeas = (TTree*) fMeas->Get("vol1");
  int meas_event_id;
  float meas_x;
  float meas_y;
  float meas_z;
  tMeas->SetBranchAddress("true_x",&meas_x);
  tMeas->SetBranchAddress("true_y",&meas_y);
  tMeas->SetBranchAddress("true_z",&meas_z);
  tMeas->SetBranchAddress("event_nr",&meas_event_id);

  // setup hits
  TFile* fHits = new TFile(TString(dir + "hits.root"));
  TTree* tHits = (TTree*) fHits->Get("hits");
  UInt_t hit_event_id;
  ULong64_t hit_particle_id;
  float hit_x;
  float hit_y;
  float hit_z;
  tHits->SetBranchAddress("event_id",&hit_event_id);
  tHits->SetBranchAddress("particle_id",&hit_particle_id);
  tHits->SetBranchAddress("tx",&hit_x);
  tHits->SetBranchAddress("ty",&hit_y);
  tHits->SetBranchAddress("tz",&hit_z);
  tHits->BuildIndex("event_id","int(tx*1048576)");
  
  // map indices of first measurements per event
  std::map<int,int> mapFirstMeaurementIndex;
  int previous_event_nr = -1;
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if (meas_event_id==previous_event_nr) continue;
    previous_event_nr = meas_event_id;
    mapFirstMeaurementIndex[meas_event_id] = im;
  }

  TH1D* hMcPtPi = new TH1D("hMcPtPi","",100,0.,1.);
  TH1D* hMcPtPr = new TH1D("hMcPtPr","",100,0.,1.);
  TH1D* hRcPtPi = new TH1D("hRcPtPi","",100,0.,1.);
  TH1D* hRcPtPr = new TH1D("hRcPtPr","",100,0.,1.);

  for (int ev=0;ev<tPart->GetEntries();ev++){
    tPart->GetEntry(ev);
    for (int ip = 0; ip<part_pdg->size(); ip++){
      int pdg = part_pdg->at(ip);
      float vz = part_vz->at(ip);
      float pt = part_pt->at(ip);
      float eta = part_eta->at(ip);
      if (abs(eta-etaMean)<etaDif && abs(vz)<1.) {
        if (abs(pdg)== 211) hMcPtPi->Fill(pt);
        if (abs(pdg)==2212) hMcPtPr->Fill(pt);
      }
    }
  } // event loop

  for (int is=0;is<tSeeds->GetEntries();is++){
    tSeeds->GetEntry(is);
    int shift = mapFirstMeaurementIndex[seed_event_id];
    tMeas->GetEntry(shift + meas_id_1);
    tHits->GetEntryWithIndex(seed_event_id, int(meas_x*1048576));
    if (abs(meas_y-hit_y)>1e-10) printf("WARNING Wrong hit index\n");
    ULong64_t hit_particle_id1 = hit_particle_id;
    tMeas->GetEntry(shift + meas_id_2);
    tHits->GetEntryWithIndex(seed_event_id, int(meas_x*1048576));
    if (abs(meas_y-hit_y)>1e-10) printf("WARNING Wrong hit index\n");
    ULong64_t hit_particle_id2 = hit_particle_id;
    tMeas->GetEntry(shift + meas_id_3);
    tHits->GetEntryWithIndex(seed_event_id, int(meas_x*1048576));
    if (abs(meas_y-hit_y)>1e-10) printf("WARNING Wrong hit index\n");
    ULong64_t hit_particle_id3 = hit_particle_id;

    if (hit_particle_id1!=hit_particle_id2 || hit_particle_id2!=hit_particle_id3) {
      // fake seed
      continue;
    }

    ActsFatras::Barcode barcode(hit_particle_id1);
    int ip = barcode.particle()-1;
    tPart->GetEntry(seed_event_id);
    if (ip>part_pdg->size()) {
      printf("strange particle id: event=%d part=%d\n", seed_event_id, ip);
      continue;
    }
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
  hMcPtPr->Draw();

  new TCanvas;
  hRcPtPi->Divide(hRcPtPi, hMcPtPi, 1, 1, "B");
  hRcPtPi->Draw();
  
  new TCanvas;
  hRcPtPr->Divide(hRcPtPr, hMcPtPr, 1, 1, "B");
  hRcPtPr->Draw();

  TFile* f = new TFile(TString(dir + "seed_efficiency.root"),"update");
  hMcPtPi->Write(Form("hMcPtPi%.0f",etaMean*10));
  hMcPtPr->Write(Form("hMcPtPr%.0f",etaMean*10));
  hRcPtPi->Write(Form("hEffPtPi%.0f",etaMean*10));
  hRcPtPr->Write(Form("hEffPtPr%.0f",etaMean*10));
  f->Close();

}
