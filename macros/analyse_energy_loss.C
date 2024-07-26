#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "ActsFatras/EventData/Barcode.hpp"
#include "analyse_trackstates.C"

//void analyse_energy_loss(TString dir = "acts_pi_16_fixedPt", double etaMean = 1.6, double etaDif = 0.1){
void analyse_energy_loss(TString dir = "acts_pr_16_fixedPt", double etaMean = 1.6, double etaDif = 0.1){
  dir.Append("/");

  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadBottomMargin(0.11);

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
  auto part_phi = new std::vector<float>;
  tPart->SetBranchAddress("event_id",&part_event_id);
  tPart->SetBranchAddress("particle_id",&part_id);
  tPart->SetBranchAddress("particle_type",&part_pdg);
  tPart->SetBranchAddress("vz",&part_vz);
  tPart->SetBranchAddress("px",&part_px);
  tPart->SetBranchAddress("py",&part_py);
  tPart->SetBranchAddress("pz",&part_pz);  
  tPart->SetBranchAddress("pt",&part_pt);  
  tPart->SetBranchAddress("eta",&part_eta);  
  tPart->SetBranchAddress("phi",&part_phi);  

  TFile* f = new TFile(dir + "trackstates.root");
  f->ls();
  TTree* t = (TTree*) f->Get("trackstates");
  SetBranchAddresses(t);

  TH1D* h = new TH1D("hDeltaPoverP","",200,0,0.4);
  for (int it=0; it<t->GetEntries(); it++){
    t->GetEntry(it);
    printf("%d %d %f\n",m_eventNr, m_trackNr, m_t_eQOP->at(0));
    double qop = m_t_eQOP->at(0);
    double p = 1/qop;
    double pMC = 0.9;
    h->Fill(-(p-pMC)/p);
  }

  
  new TCanvas;
  h->SetTitle(";#Delta p/p;");
  h->SetLineWidth(2);
  h->SetLineColor(kBlue);

  h->SetLabelSize(0.045,"XY");
  h->SetTitleSize(0.045,"XY");
  h->SetTitleOffset(1.1,"X");
  h->Draw();
  gPad->Print("delta_p_over_p_for_pr.png");
}
