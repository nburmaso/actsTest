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

void analyse_measurements(){
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

  TH1D* hLoc0 = new TH1D("hLoc0","",1000,-10,10);
  TH1D* hSurfaceID0 = new TH1D("hSurfaceID0","",1000,0,1000);
  TH1D* hSurfaceID4 = new TH1D("hSurfaceID4","",1000,0,1000);

  TH1D* hSurfaceIDLayer1 = new TH1D("hSurfaceID0Layer1","",1000,0,1000);
  TH1D* hSurfaceIDLayer2 = new TH1D("hSurfaceID0Layer2","",1000,0,1000);

  TH2D* hL1vsL4 = new TH2D("hL1vsL4","",1000,0,1000,1000,0,1000);
  TH1D* hL4minusL0 = new TH1D("hL4minusL0","",10,-10,10);
  TH1D* hL5minusL1 = new TH1D("hL5minusL1","",10,-10,10);
  TH1D* hL6minusL2 = new TH1D("hL6minusL2","",10,-10,10);

  TH1D* hL1minusL0 = new TH1D("hL1minusL0","",10,-10,10);
  TH1D* hL2minusL1 = new TH1D("hL2minusL1","",10,-10,10);
  TH1D* hL2minusL0 = new TH1D("hL2minusL0","",10,-10,10);
  TH1D* hL5minusL0 = new TH1D("hL5minusL0","",10,-10,10);
  TH1D* hL6minusL0 = new TH1D("hL6minusL0","",10,-10,10);

  int prev_event_id = -1;
  int l0 = -1;
  int l1 = -1;
  int l2 = -1;
  int l4 = -1;
  int l5 = -1;
  int l6 = -1;
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if (m_event_id!=prev_event_id) {
      l0=-1;
      l1=-1;
      l2=-1;
      l4=-1;
      l5=-1;
      l6=-1;
      prev_event_id = m_event_id;
    }
    if (m_layer_id==2) continue;
    if (m_layer_id==38) continue;    
    if ((m_layer_id-3)%7==3) continue;
    // printf("%f\n",m_true_loc0);
    hLoc0->Fill(m_true_loc0);
    if ((m_layer_id-3)/7==0) hSurfaceID0->Fill(m_surface_id);
    if ((m_layer_id-3)/7==4) hSurfaceID4->Fill(m_surface_id);
    if ((m_layer_id-3)==0) {
      hSurfaceIDLayer1->Fill(m_surface_id);
      l0=m_surface_id;
    }
    if ((m_layer_id-3)==1) {
      l1=m_surface_id;
      if (l0!=-1) hL1minusL0->Fill(m_surface_id-l0);
    }
    if ((m_layer_id-3)==2) {
      l2=m_surface_id;
      if (l0!=-1) hL2minusL0->Fill(m_surface_id-l0);
      if (l1!=-1) hL2minusL1->Fill(m_surface_id-l1);
    }
    if ((m_layer_id-3)==4) {
      l5=m_surface_id;
    }
    if ((m_layer_id-3)==5) {
      l5=m_surface_id;
      if (l0!=-1) hL5minusL0->Fill(m_surface_id-l0);
    }
    if ((m_layer_id-3)==6) {
      l6=m_surface_id;
      if (l0!=-1) hL6minusL0->Fill(m_surface_id-l0);
    }

    if ((m_layer_id-3)==4) {
      hSurfaceIDLayer2->Fill(m_surface_id);
      printf("%d\n",l0);
      if (l0!=-1) hL4minusL0->Fill(m_surface_id-l0);
    }
    if ((m_layer_id-3)==5) {
      if (l1!=-1) hL5minusL1->Fill(m_surface_id-l1);
    }
    if ((m_layer_id-3)==6) {
      if (l2!=-1) hL6minusL2->Fill(m_surface_id-l2);
    }
  }
  new TCanvas;
  hLoc0->Draw();
  new TCanvas;
  hSurfaceID0->Draw();
  new TCanvas;
  hSurfaceID4->Draw();
  new TCanvas;
  hSurfaceIDLayer1->Draw();
  new TCanvas;
  hSurfaceIDLayer2->Draw();
  new TCanvas;
  hL1vsL4->Draw();
  new TCanvas;
  hL4minusL0->Draw();
  new TCanvas;
  hL5minusL1->Draw();
  new TCanvas;
  hL6minusL2->Draw();
  new TCanvas;
  hL1minusL0->Draw();
  new TCanvas;
  hL2minusL0->Draw();
  new TCanvas;
  hL5minusL0->Draw();
  new TCanvas;
  hL6minusL0->Draw();
  new TCanvas;
  hL2minusL1->Draw();
}

