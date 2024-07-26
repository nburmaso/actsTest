#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "style.h"

void draw_pulls(TString dir = "notpc_pi_16", TString pid = "Pi", double etaMean = 1.6, double minPtMC=0.25, double maxPtMC=0.3){
  dir.Append("/");
//  gStyle->SetStatFontSize(0.08);
  gStyle->SetStatH(0.15);
  gStyle->SetStatW(0.2);
  gStyle->SetStatFormat("6.3g");
  gStyle->SetOptFit(101);

  TFile* f = new TFile(dir + "tracking_efficiency.root");
  TH2D* hResQOPvsPt  = (TH2D*) f->Get(Form("hResQOPvsPt%s%.0f",pid.Data(),etaMean*10));
  TH2D* hResDvsPt    = (TH2D*) f->Get(Form("hResDvsPt%s%.0f",pid.Data(),etaMean*10));
  TH2D* hResZvsPt    = (TH2D*) f->Get(Form("hResZvsPt%s%.0f",pid.Data(),etaMean*10));
  TH2D* hPullQOPvsPt = (TH2D*) f->Get(Form("hPullQOPvsPt%s%.0f",pid.Data(),etaMean*10));
  TH2D* hPullDvsPt   = (TH2D*) f->Get(Form("hPullDvsPt%s%.0f",pid.Data(),etaMean*10));
  TH2D* hPullZvsPt   = (TH2D*) f->Get(Form("hPullZvsPt%s%.0f",pid.Data(),etaMean*10));

  int iMinBin = hResQOPvsPt->GetXaxis()->FindFixBin(minPtMC+0.0001);
  int iMaxBin = hResQOPvsPt->GetXaxis()->FindFixBin(maxPtMC-0.0001);

  TCanvas* c1 = new TCanvas("c1","c1",1800,800);
  c1->Divide(3,2,0.001,0.001);
  
  TF1* fGaus = new TF1("fGaus","gaus", -1, 1);
  TH1D* hResQOP  = hResQOPvsPt->ProjectionY("hResQOP",iMinBin,iMaxBin);
  TH1D* hResD    = hResDvsPt->ProjectionY("hResD",iMinBin,iMaxBin);
  TH1D* hResZ    = hResZvsPt->ProjectionY("hResZ",iMinBin,iMaxBin);
  TH1D* hPullQOP = hPullQOPvsPt->ProjectionY("hPullQOP",iMinBin,iMaxBin);
  TH1D* hPullD   = hPullDvsPt->ProjectionY("hPullD",iMinBin,iMaxBin);
  TH1D* hPullZ   = hPullZvsPt->ProjectionY("hPullZ",iMinBin,iMaxBin);

  c1->cd(1);
  SetPad(gPad);
  SetHisto(hResQOP,Form(";(q/p)^{Rec} - (q/p)^{MC} (1/GeV)"));
  hResQOP->Draw();
  hResQOP->Fit(fGaus,"LQ","");
  c1->cd(2);
  SetPad(gPad);
  SetHisto(hResD,Form(";d (cm)"));
  hResD->Draw();
  hResD->Fit(fGaus,"LQ","");
  c1->cd(3);
  SetPad(gPad);
  SetHisto(hResZ,Form(";z (cm)"));
  hResZ->Draw();
  hResZ->Fit(fGaus,"LQ","");
  c1->cd(4);
  SetPad(gPad);
  SetHisto(hPullQOP,Form(";pull q/p"));
  hPullQOP->Draw();
  hPullQOP->Fit(fGaus,"LQ","",-4,4);
  c1->cd(5);
  SetPad(gPad);
  SetHisto(hPullD,Form(";pull d"));
  hPullD->Draw();
  hPullD->Fit(fGaus,"LQ","",-4,4);
  c1->cd(6);
  SetPad(gPad);
  SetHisto(hPullZ,Form(";pull z"));
  hPullZ->Draw();
  hPullZ->Fit(fGaus,"LQ","",-4,4);
  c1->Print(dir+Form("pulls_%s_%.0f_%.0f.png", pid.Data(), etaMean*10, minPtMC*1000));
}
