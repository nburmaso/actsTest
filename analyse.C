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
  TFile* f = new TFile("trackstates_ckf.root");
  f->ls();
  vector<float> t_x;     vector<float>* ptr_t_x = &t_x;

  UInt_t track_nr;
  TTree* t = (TTree*) f->Get("trackstates");
  t->Print();
  t->SetBranchAddress("t_x",&ptr_t_x);
  t->SetBranchAddress("track_nr",&track_nr);
  for (int itr = 0;itr<t->GetEntries();itr++){
    t->GetEntry(itr);
    printf("%d\n",track_nr);
    for (UInt_t i=0;i<t_x.size();i++){
      printf("  %f\n",t_x[i]);
    }
  }
}