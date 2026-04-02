#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
//#include "MpdMCTrack.h"
//#include "MpdFwdPoint.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"


const int nStations = 5;
const int nLayersPerStation = 9;
const int shift = 2;

void analyse_measurements(int selected_station = 4){
  gStyle->SetOptStat(0);
  TFile* fMeas = new TFile("../build/test/measurements.root");
  TTree* tMeas = (TTree*) fMeas->Get("measurements");
  // tMeas->Print();
  int32_t m_event_id;
  int32_t m_volume_id;
  int32_t m_layer_id;
  int32_t m_surface_id;
  float m_true_loc0;
  tMeas->SetBranchAddress("event_nr",&m_event_id);
  tMeas->SetBranchAddress("volume_id",&m_volume_id);
  tMeas->SetBranchAddress("layer_id",&m_layer_id);
  tMeas->SetBranchAddress("surface_id",&m_surface_id);
  tMeas->SetBranchAddress("true_loc0",&m_true_loc0);

  TH1D* hStrawDiff[nLayersPerStation];
  for (int i=0;i<nLayersPerStation;i++){
    hStrawDiff[i] = new TH1D(Form("hStrawDiff%d",i),Form("Layer%d - Layer0",i),20,-10,10);
  }

  int prev_event_id = -1;
  int straw0 = -1;
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if (m_event_id!=prev_event_id) {
      straw0 = -1;
      prev_event_id = m_event_id;
    }
    int station = (m_layer_id-shift)/nLayersPerStation;
    int layer = (m_layer_id-shift)%nLayersPerStation;
    int straw = m_surface_id;
    if (station!=selected_station) continue;
    if (layer==0) straw0 = straw;
    if (straw0<0) continue;
    hStrawDiff[layer]->Fill(straw-straw0);
  }

  TCanvas* cStrawDiff = new TCanvas(Form("straw_diff_station%d",selected_station),Form("straw_diff_station%d",selected_station),1900,1000);
  cStrawDiff->Divide(3,3,0.001,0.001);
  for (int i=0;i<nLayersPerStation;i++){
    cStrawDiff->cd(i+1);
    hStrawDiff[i]->Draw();
  }
  cStrawDiff->Print(".png");
}

