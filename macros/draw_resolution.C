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


void draw_resolution(TString pid = "Pi", TString dir1 = "notpc_pi_16", TString dir2 = "notpc_pi_19", TString dir3 = "notpc_pi_22"){
  dir1.Append("/");
  dir2.Append("/");
  dir3.Append("/");

  TFile* f1 = new TFile(dir1 + "resolution.root");
  TFile* f2 = new TFile(dir2 + "resolution.root");
  TFile* f3 = new TFile(dir3 + "resolution.root");
  TGraph* gRes16 = (TGraph*) f1->Get(Form("gRes%s16",pid.Data()));
  TGraph* gRes19 = (TGraph*) f2->Get(Form("gRes%s19",pid.Data()));
  TGraph* gRes22 = (TGraph*) f3->Get(Form("gRes%s22",pid.Data()));

  const int nPtBins = 25;
  double pt[nPtBins]={0.005,0.007,0.01,0.015,0.02,0.03,0.04,0.05,0.07,0.1,0.12,0.15,0.2,0.22, 0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,5}; // GeV
  double delta_pt_ms[nPtBins];
  double delta_pt_res[nPtBins];
  double delta_pt[nPtBins];

  double delta2_pt_ms[nPtBins];
  double delta2_pt_res[nPtBins];
  double delta2_pt[nPtBins];
  
  const int nEta = 3;
  double veta[nEta]={1.6, 1.9, 2.2};
  double theta[nEta];
  for (int i=0;i<nEta;i++) theta[i]= 2*atan(exp(-veta[i]));

  TGraph* vGraphRes[nEta];

  // double m = 0.10566;  // GeV
  // double m = 0.13957;  // GeV
  // double m = 0.938;    // GeV
  double m = pid.Contains("Pi") ? 0.13957 : 0.938;
  double B = 0.5;             // T
  double sigma_rf = 0.005;     // cm
  double d_over_X0 = 0.001067; // 100 us of silicon
  double L0 = 0.9; // m
  double N= 4;
  fMS = new TF1("fMS","[0]/sqrt(([0]+1)*([0]-1))*0.0136/(0.3*x/sqrt(x*x+[3]*[3])*[1])*sqrt(([0]+1)*[2])*(1+0.038*log([2]))",0.001,1.1);
  fRes = new TF1("fRes","[0]*x",0.001,1.1);
  TF1* fSum = new TF1("fSum",quadratic_sum,0.001,1.1,5);
  fMS->SetLineWidth(1);
  fRes->SetLineWidth(1);
  fSum->SetLineWidth(1);
  // for (int ieta = 0; ieta < nEta; ieta++){
  //   double eta = veta[ieta];
  //   for (int i = 0;i<nPtBins;i++){
  //     double N = 4;
  //     double mat = (N+1)*d_over_X0;
  //     double tan_half_theta = TMath::Exp(-eta);
  //     double theta = 2 * TMath::ATan(tan_half_theta);
  //     double cost = cos(theta);
  //     double p = pt[i]/sin(theta);
  //     double beta = p/sqrt(p*p+m*m);
  //     double L = L0*tan(theta); // m
  //     delta_pt_ms[i]=(N/sqrt((N+1)*(N-1)))*(0.0136/(0.3*beta*B*L))*sqrt(mat/cost)*(1+0.038*log(d_over_X0/cost));
  //     delta_pt_res[i]=sigma_rf*pt[i]/(0.3*B*L*L)*sqrt(720*N*N*N/((N-1)*(N+1)*(N+2)*(N+3)))/100;
  //     delta_pt[i]= sqrt(delta2_pt_ms[i]*delta2_pt_ms[i]+delta2_pt_res[i]*delta2_pt_res[i]);
  //     printf("%f\n",delta_pt_ms[i]);
  //   }
    
  //   vGraphRes[ieta] = new TGraph(nPtBins,pt,delta_pt_res);
  //   vGraphRes[ieta]->SetMarkerStyle(kFullCircle);
  //   vGraphRes[ieta]->SetLineWidth(2);
  //   vGraphRes[ieta]->SetLineStyle(9);
  // }
  // vGraphRes[0]->SetLineColor(kBlue);
  // vGraphRes[1]->SetLineColor(kMagenta);
  // vGraphRes[2]->SetLineColor(kRed);

  // vGraphRes[0]->Draw();
  // return;

  TCanvas* c = new TCanvas("c", "c", 1200, 900);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.02);
  TH1F* hFrame = gPad->DrawFrame(0.,0.,1.1,0.15);
  hFrame->SetTitle(";p_{T}^{MC} (GeV/c); Resolution");
  hFrame->GetXaxis()->SetTitleOffset(1.2);

  fMS->SetParameters(N, B*L0*tan(theta[0]), d_over_X0/cos(theta[0]), m*sin(theta[0]));
  fMS->SetLineColor(kBlue-10);
  fMS->DrawCopy("same");
  fRes->SetParameter(0,sigma_rf/(0.3*B*L0*L0*tan(theta[0])*tan(theta[0]))*sqrt(720*N*N*N/((N-1)*(N+1)*(N+2)*(N+3)))/100);
  fRes->SetLineColor(kBlue-10);
  fRes->DrawCopy("same");
  fSum->SetParameters(fRes->GetParameter(0),fMS->GetParameter(0),fMS->GetParameter(1),fMS->GetParameter(2),fMS->GetParameter(3));
  fSum->SetLineColor(kBlue);
  fSum->DrawCopy("same");

  fMS->SetParameters(N, B*L0*tan(theta[1]), d_over_X0/cos(theta[1]), m*sin(theta[1]));
  fMS->SetLineColor(kMagenta-10);
  fMS->DrawCopy("same");
  fRes->SetParameter(0,sigma_rf/(0.3*B*L0*L0*tan(theta[1])*tan(theta[1]))*sqrt(720*N*N*N/((N-1)*(N+1)*(N+2)*(N+3)))/100);
  fRes->SetLineColor(kMagenta-10);
  fRes->DrawCopy("same");
  fSum->SetParameters(fRes->GetParameter(0),fMS->GetParameter(0),fMS->GetParameter(1),fMS->GetParameter(2),fMS->GetParameter(3));
  fSum->SetLineColor(kMagenta);
  fSum->DrawCopy("same");
  
  fMS->SetParameters(N, B*L0*tan(theta[2]), d_over_X0/cos(theta[2]), m*sin(theta[2]));
  fMS->SetLineColor(kRed-10);
  fMS->DrawCopy("same");
  fRes->SetParameter(0,sigma_rf/(0.3*B*L0*L0*tan(theta[2])*tan(theta[2]))*sqrt(720*N*N*N/((N-1)*(N+1)*(N+2)*(N+3)))/100);
  fRes->SetLineColor(kRed-10);
  fRes->DrawCopy("same");
  fSum->SetParameters(fRes->GetParameter(0),fMS->GetParameter(0),fMS->GetParameter(1),fMS->GetParameter(2),fMS->GetParameter(3));
  fSum->SetLineColor(kRed);
  fSum->DrawCopy("same");

  // vGraphRes[0]->Draw();
  // vGraphRes[1]->Draw("same");
  // vGraphRes[2]->Draw("same");

  gRes16->SetLineWidth(2); gRes16->SetLineColor(kBlue);
  gRes19->SetLineWidth(2); gRes19->SetLineColor(kMagenta);
  gRes22->SetLineWidth(2); gRes22->SetLineColor(kRed);

  if (dir1.Contains("acts") && pid.Contains("Pr")) {
    gRes16->RemovePoint(0);
    gRes19->RemovePoint(0);
    gRes22->RemovePoint(0);
  }

  gRes16->Draw("same");
  gRes19->Draw("same");
  gRes22->Draw("same");

  TLegend* legend = new TLegend(0.45,0.58,0.65,0.75);
  legend->SetBorderSize(0);
  legend->AddEntry(gRes22,"#eta = 2.2","l");
  legend->AddEntry(gRes19,"#eta = 1.9","l");
  legend->AddEntry(gRes16,"#eta = 1.6","l");
  legend->Draw();
  
  gPad->Print(dir1 + Form("resolution%s.png",pid.Data()));
}

