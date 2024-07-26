#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "map"
#include "vector"
#include "TEllipse.h"
#include "TGraphErrors.h"

using namespace std;

void analyse(){
  TFile* f = new TFile("build/hits.root");
  float tx;
  float ty;
  float tz;

  UInt_t track_nr;
  TTree* t = (TTree*) f->Get("hits");
  t->SetBranchAddress("tx",&tx);
  t->SetBranchAddress("ty",&ty);
  t->SetBranchAddress("tz",&tz);

  double minX = -155;
  double maxX =  155;
  double minY = -155;
  double maxY =  155;
  int nX = (maxX-minX)*10;
  int nY = (maxY-minY)*10;

  TH2D* hXY = new TH2D("hXY","",nX,minX,maxX,nY,minY,maxY);
  for (int i = 0;i<t->GetEntries();i++){
    t->GetEntry(i);
    if (tz>2000 && tz<2150){
    //if (tz>2900 && tz<3100){
    //if (tz>1500 && tz<1660){
    //if (tz>1660 && tz<1720){
      hXY->Fill(tx/10.,ty/10.);
    }
  }

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("c1","c1",950,950);
  gPad->SetRightMargin(0.12);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.10);
  hXY->SetTitle("TPC flange implementated in ACTS;cm;cm");
  hXY->Draw("colz");
  gPad->Print("flange.png");
}