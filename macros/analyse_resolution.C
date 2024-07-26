#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "style.h"

void analyse_resolution(TString dir = "notpc_pi_16_200_5", TString pid = "Pi", double etaMean = 1.6){
  dir.Append("/");
//  gStyle->SetStatFontSize(0.08);
  gStyle->SetStatH(0.15);
  gStyle->SetStatW(0.2);
  gStyle->SetStatFormat("6.3g");
  gStyle->SetOptFit(101);

  TFile* f = new TFile(dir + "tracking_efficiency.root");
  TH2D* hPtResVsPtMC = (TH2D*) f->Get(Form("hPtResVsPt%s%.0f",pid.Data(), etaMean*10));

  int nbins = hPtResVsPtMC->GetNbinsX();
  
  TCanvas* c1 = new TCanvas("c1","c1",1800,800);
  c1->Divide(4,2,0.001,0.02);
  
  
  const int nPtBins = 18;
  double vPt[nPtBins];
  double vRes[nPtBins];

  TF1* fGaus = new TF1("fGaus","gaus", -1, 1);
  int imin = 3;
  for (int ibin = imin; ibin<=nbins; ibin++){
    double minPtMC = 1000*hPtResVsPtMC->GetXaxis()->GetBinLowEdge(ibin);
    double maxPtMC = 1000*hPtResVsPtMC->GetXaxis()->GetBinUpEdge(ibin);
    TH1D* hProj = hPtResVsPtMC->ProjectionY(Form("hRes_%.0f_%.0f",minPtMC,maxPtMC),ibin,ibin);
    hProj->Fit(fGaus,"LQ0","",-0.2,0.2);
    vPt[ibin-imin] = hPtResVsPtMC->GetXaxis()->GetBinCenter(ibin);
    vRes[ibin-imin] = fGaus->GetParameter(2);
    //vRes[ibin-imin] = hProj->GetRMS();
    if (ibin%2==1) continue;
    c1->cd(ibin/2-1);
    SetPad(gPad);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.15);
    SetHisto(hProj,Form("p_{T}: %.0f - %.0f MeV;p_{T} resolution: (p_{T}^{Rec} - p_{T}^{MC})/p_{T}^{MC}",minPtMC,maxPtMC));
    hProj->SetTitleOffset(1.3);
    hProj->GetXaxis()->SetRangeUser(-0.7,0.7);
    hProj->Draw();
    hProj->GetListOfFunctions()->At(0)->Draw("same");
    
  }
  c1->Print(dir + "res.png");

  for (int i=0;i<nPtBins;i++){
    printf("%d %f %f\n",i,vPt[i],vRes[i]);
  }
  TGraph* g = new TGraph(nPtBins,vPt,vRes);
  TFile* fg = new TFile(dir + "resolution.root","update");
  g->Write(Form("gRes%s%.0f",pid.Data(), etaMean*10));
  fg->Close();
}
