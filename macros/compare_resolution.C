#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "MpdMCTrack.h"
#include "MpdFwdPoint.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"

TF1* fMS;
TF1* fRes;
double quadratic_sum(double* x, double* par){
   fRes->SetParameter(0,par[0]);
   fMS->SetParameter(0,par[1]);
   fMS->SetParameter(1,par[2]);
   fMS->SetParameter(2,par[3]);
   fMS->SetParameter(3,par[4]);

   return sqrt(pow(fRes->Eval(x[0]),2) + pow(fMS->Eval(x[0]),2.));
}


void compare_resolution(TString pid = "Pr"){
  bool isPi = pid.Contains("Pi");
  TFile* f1 = new TFile(Form("notpc_%s_16/resolution.root",isPi?"pi":"pr"));
  TFile* f2 = new TFile(Form("notpc_%s_19/resolution.root",isPi?"pi":"pr"));
  TFile* f3 = new TFile(Form("notpc_%s_22/resolution.root",isPi?"pi":"pr"));
  TGraph* gRes16 = (TGraph*) f1->Get(Form("gRes%s16",pid.Data()));
  TGraph* gRes19 = (TGraph*) f2->Get(Form("gRes%s19",pid.Data()));
  TGraph* gRes22 = (TGraph*) f3->Get(Form("gRes%s22",pid.Data()));

  TFile* fRoc1 = new TFile(Form("noframe_%s_16/resolution.root",isPi?"pi":"pr"));
  TFile* fRoc2 = new TFile(Form("noframe_%s_19/resolution.root",isPi?"pi":"pr"));
  TFile* fRoc3 = new TFile(Form("noframe_%s_22/resolution.root",isPi?"pi":"pr"));
  TGraph* gResRoc16 = (TGraph*) fRoc1->Get(Form("gRes%s16",pid.Data()));
  TGraph* gResRoc19 = (TGraph*) fRoc2->Get(Form("gRes%s19",pid.Data()));
  TGraph* gResRoc22 = (TGraph*) fRoc3->Get(Form("gRes%s22",pid.Data()));

  TFile* fFull1 = new TFile(Form("acts_%s_16/resolution.root",isPi?"pi":"pr"));
  TFile* fFull2 = new TFile(Form("acts_%s_19/resolution.root",isPi?"pi":"pr"));
  TFile* fFull3 = new TFile(Form("acts_%s_22/resolution.root",isPi?"pi":"pr"));
  TGraph* gResFull16 = (TGraph*) fFull1->Get(Form("gRes%s16",pid.Data()));
  TGraph* gResFull19 = (TGraph*) fFull2->Get(Form("gRes%s19",pid.Data()));
  TGraph* gResFull22 = (TGraph*) fFull3->Get(Form("gRes%s22",pid.Data()));


  TCanvas* c = new TCanvas("c", "c", 1200, 900);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.02);
  TH1F* hFrame = gPad->DrawFrame(0.,0.,1.1,0.15);
  hFrame->SetTitle(";p_{T}^{MC} (GeV/c); Resolution");
  hFrame->GetXaxis()->SetTitleOffset(1.2);

  gRes16->SetLineWidth(2); gRes16->SetLineColor(kBlue-10);
  gRes19->SetLineWidth(2); gRes19->SetLineColor(kMagenta-10);
  gRes22->SetLineWidth(2); gRes22->SetLineColor(kRed-10);

  gResRoc16->SetLineWidth(2); gResRoc16->SetLineColor(kBlue-9);
  gResRoc19->SetLineWidth(2); gResRoc19->SetLineColor(kMagenta-9);
  gResRoc22->SetLineWidth(2); gResRoc22->SetLineColor(kRed-9);

  gResFull16->SetLineWidth(2); gResFull16->SetLineColor(kBlue);
  gResFull19->SetLineWidth(2); gResFull19->SetLineColor(kMagenta);
  gResFull22->SetLineWidth(2); gResFull22->SetLineColor(kRed);

  gRes16->Draw("same");
  gRes19->Draw("same");
  gRes22->Draw("same");

  gResRoc16->RemovePoint(0);
  // gResRoc19->RemovePoint(0);
  // gResRoc22->RemovePoint(0);

  gResRoc16->Draw("same");
  gResRoc19->Draw("same");
  gResRoc22->Draw("same");

  gResFull16->RemovePoint(0);
  gResFull19->RemovePoint(0);
  gResFull22->RemovePoint(0);

  gResFull16->RemovePoint(0);
  gResFull19->RemovePoint(0);
  gResFull22->RemovePoint(0);

  // gResFull16->Draw("same");
  // gResFull19->Draw("same");
  // gResFull22->Draw("same");

  TLegend* legend = new TLegend(0.35,0.68,0.60,0.85);
  legend->SetBorderSize(0);
  legend->AddEntry(gResFull22,"#eta = 2.2","l");
  legend->AddEntry(gResFull19,"#eta = 1.9","l");
  legend->AddEntry(gResFull16,"#eta = 1.6","l");
  legend->Draw();
  
  gPad->Print(Form("resolution%s.png",pid.Data()));
}

