#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "ActsFatras/EventData/Barcode.hpp"
#include "analyse_summary.C"

void analyse_tracking_efficiency(TString dir = "acts_pi_16", double etaMean = 1.6, double etaDif = 0.1){
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

  TFile* fTrack = new TFile(TString(dir + "tracksummary.root"));
  TTree* tTrack = (TTree*) fTrack->Get("tracksummary");
  SetBranchAddresses(tTrack);
  int nEvents = tTrack->GetEntries();

  TH1D* hMcPtPi = new TH1D("hMcPtPi","",100,0.,1.);
  TH1D* hMcPtPr = new TH1D("hMcPtPr","",100,0.,1.);
  TH1D* hRcPtPi = new TH1D("hRcPtPi","",100,0.,1.);
  TH1D* hRcPtPr = new TH1D("hRcPtPr","",100,0.,1.);

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
    } // particle loop
  } // event loop
  
  TVector3 v;
  TVector3 vMC;
  for (int ev=0; ev<tTrack->GetEntries(); ev++){
    tTrack->GetEntry(ev);
    tPart->GetEntry(ev);

    for (int it=0; it<m_majorityParticleId->size(); it++){
      if (!m_hasFittedParams->at(it)) continue;
      double qp = m_eQOP_fit->at(it);
      double theta = m_eTHETA_fit->at(it);
      double phi = m_ePHI_fit->at(it);
      v.SetMagThetaPhi(1./qp, theta, phi);
      double ptRC = v.Pt();

      ActsFatras::Barcode barcode(m_majorityParticleId->at(it));
      int ip = barcode.particle()-1;
      if (ip==65534) continue;
      int pdg = part_pdg->at(ip);
      float vz = part_vz->at(ip);
      float ptMC = part_pt->at(ip);
      float etaMC = part_eta->at(ip);
      float phiMC = part_phi->at(ip);
      vMC.SetPtEtaPhi(ptMC, etaMC, phiMC);
      float qpMC = 1/vMC.Mag();
      if (abs(etaMC-etaMean)<etaDif && abs(vz)<1.) {
        if (abs(pdg)== 211) {
          hRcPtPi->Fill(ptMC);
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
  hMcPtPi->Draw();

  new TCanvas;
  hMcPtPr->Draw();

  new TCanvas;
  hRcPtPi->Divide(hRcPtPi, hMcPtPi, 1, 1, "B");
  hRcPtPi->Draw();
  
  new TCanvas;
  hRcPtPr->Divide(hRcPtPr, hMcPtPr, 1, 1, "B");
  hRcPtPr->Draw();

  new TCanvas;
  hPtResVsPtPi->Draw();

  new TCanvas;
  hPtResVsPtPr->Draw();

  TFile* f = new TFile(TString(dir + "tracking_efficiency.root"),"update");
  hMcPtPi->Write(Form("hMcPtPi%.0f",etaMean*10));
  hMcPtPr->Write(Form("hMcPtPr%.0f",etaMean*10));
  hRcPtPi->Write(Form("hEffPtPi%.0f",etaMean*10));
  hRcPtPr->Write(Form("hEffPtPr%.0f",etaMean*10));
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
