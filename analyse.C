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
  TFile* f = new TFile("hits.root");
  float tx;
  float ty;
  float tz;

  UInt_t track_nr;
  TTree* t = (TTree*) f->Get("hits");
  t->SetBranchAddress("tx",&tx);
  t->SetBranchAddress("ty",&ty);
  t->SetBranchAddress("tz",&tz);

  TH2D* hXY = new TH2D("hXY","",1000,-1500,1500,1000,-1500,1500);
  for (int i = 0;i<t->GetEntries();i++){
    t->GetEntry(i);
    // if (tz>2000 && tz<2150){
    //if (tz>1500 && tz<1660){
    if (tz>1660 && tz<1720){
      hXY->Fill(tx,ty);
    }
  }

  new TCanvas;
  hXY->Draw("colz");
}