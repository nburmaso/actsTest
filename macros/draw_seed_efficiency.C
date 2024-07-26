#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"

void draw_seed_efficiency(TString pid = "pi", TString dir1 = "acts_pi_16", TString dir2 = "acts_pi_19", TString dir3 = "acts_pi_22"){
  dir1.Append("/");
  dir2.Append("/");  
  dir3.Append("/");
  bool isPi = pid.Contains("pi");

  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadBottomMargin(0.11);

  TFile* f1 = new TFile(dir1 + "seed_efficiency.root");
  TFile* f2 = new TFile(dir2 + "seed_efficiency.root");
  TFile* f3 = new TFile(dir3 + "seed_efficiency.root");
  TH1D* hEffPt1 = (TH1D*) f1->Get(isPi ? "hEffPtPi16" : "hEffPtPr16");
  TH1D* hEffPt2 = (TH1D*) f2->Get(isPi ? "hEffPtPi19" : "hEffPtPr19");
  TH1D* hEffPt3 = (TH1D*) f3->Get(isPi ? "hEffPtPi22" : "hEffPtPr22");

  hEffPt1->SetTitle(";p_{T} (GeV);Efficiency");
  hEffPt1->SetLineWidth(2);
  hEffPt2->SetLineWidth(2);
  hEffPt3->SetLineWidth(2);
  hEffPt1->SetLineColor(kBlue);
  hEffPt2->SetLineColor(kMagenta);
  hEffPt3->SetLineColor(kRed);
  
  new TCanvas;
  hEffPt1->GetYaxis()->SetRangeUser(0,1.1);
  hEffPt1->SetLabelSize(0.045,"XY");
  hEffPt1->SetTitleSize(0.045,"XY");
  hEffPt1->SetTitleOffset(1.1,"X");
  hEffPt1->Draw();
  hEffPt2->Draw("same");
  hEffPt3->Draw("same");
}
