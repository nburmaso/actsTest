#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"

void draw_mc_spectra(TString dir="mpd"){
  dir.Append("/");
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadBottomMargin(0.11);


  TFile* f1 = new TFile(dir + "tracking_efficiency.root");
  TH1D* hPt1 = (TH1D*) f1->Get("hMcPtPi22");
  TH1D* hPt2 = (TH1D*) f1->Get("hMcPtPi19");
  TH1D* hPt3 = (TH1D*) f1->Get("hMcPtPi16");
  TH1D* hPt4 = (TH1D*) f1->Get("hMcPtPr22");
  TH1D* hPt5 = (TH1D*) f1->Get("hMcPtPr19");
  TH1D* hPt6 = (TH1D*) f1->Get("hMcPtPr16");


  new TCanvas;
  hPt2->SetLabelSize(0.045,"XY");
  hPt2->SetTitleSize(0.045,"XY");
  hPt2->SetTitleOffset(1.1,"X");
  hPt2->SetLineColor(kBlue);
  hPt2->SetLineWidth(2);
  hPt2->SetTitle(";p_{T} (GeV);");
  hPt2->Draw();
  // hPt1->Draw("same");
  // hPt3->Draw("same");

  gPad->Print("ptMC_pi.png");

  new TCanvas;
  hPt5->SetLabelSize(0.045,"XY");
  hPt5->SetTitleSize(0.045,"XY");
  hPt5->SetTitleOffset(1.1,"X");
  hPt5->SetLineColor(kBlue);
  hPt5->SetLineWidth(2);
  hPt5->SetTitle(";p_{T} (GeV);");
  hPt5->Draw();
  // hPt4->Draw("same");
  // hPt6->Draw("same");

  gPad->Print("ptMC_pr.png");
}
