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

void draw(){
  TFile* f = new TFile("hits.root");
  TTree* t = (TTree*) f->Get("hits");
  float tz = 0;
  float tx = 0;
  float ty = 0;
  UInt_t event_id = 0;
  ULong64_t particle_id = 0;
  t->SetBranchAddress("event_id",&event_id);
  t->SetBranchAddress("particle_id",&particle_id);
  t->SetBranchAddress("tx",&tx);
  t->SetBranchAddress("ty",&ty);
  t->SetBranchAddress("tz",&tz);
  TH1D* hZ = new TH1D("hZ","",100,2000,3100);
  const int nMax = 1000;
  double x[nMax];
  double y[nMax];
  int n=0;
  printf("%lld\n",t->GetEntries());

  int selectedEvent = 55;

  std::map<ULong64_t,TGraph*> mapGraphHitsXY;
  for (int entry=0;entry<t->GetEntries();entry++){
    t->GetEntry(entry);
    //printf("%d",event_id);
    if (event_id!=selectedEvent) continue;
    hZ->Fill(tz);
    //printf("%llu\n",particle_id);
    if (auto graph = mapGraphHitsXY.find(particle_id); graph == mapGraphHitsXY.end()) mapGraphHitsXY[particle_id] = new TGraph();
    mapGraphHitsXY[particle_id]->AddPoint(tx,ty);
  }
  new TCanvas;
  hZ->Draw();
  new TCanvas("chits","chits",1200,1200);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.01);

  TH1F* hFrame = gPad->DrawFrame(-1300,-1300,1300,1300);
  TEllipse* elDot = new TEllipse(0,0, 20);
  TEllipse* elMin = new TEllipse(0,0, 357);
  TEllipse* elMax = new TEllipse(0,0,1300);
  elDot->SetFillStyle(0);
  elMin->SetFillStyle(0);
  elMax->SetFillStyle(0);
  elMax->Draw();
  elMin->Draw();
  elDot->Draw();

  for (auto& graph : mapGraphHitsXY){
    graph.second->SetMarkerStyle(kFullCircle);
    graph.second->SetMarkerSize(0.5);
    graph.second->Draw("lp");
  }


  ifstream fmeasurements(Form("measurements/event%09d-measurements.csv",selectedEvent));
  vector<double> xMes;
  vector<double> yMes;
  vector<double> xErr;
  vector<double> yErr;
  for( std::string line; getline(fmeasurements, line); ){
    TObjArray* objArray = TString(line).Tokenize(",");
    if (TString(objArray->At(0)->GetName()).Contains("measurement_id")) continue;
    float x = TString(objArray->At( 3)->GetName()).Atof();
    float y = TString(objArray->At( 4)->GetName()).Atof();
    float xVar = TString(objArray->At( 8)->GetName()).Atof();
    float yVar = TString(objArray->At( 9)->GetName()).Atof();
    // printf("%f %f %f %f\n",x,y,xVar,yVar);
    xMes.push_back(x);
    yMes.push_back(y);
    xErr.push_back(sqrt(xVar));
    yErr.push_back(sqrt(yVar));
    objArray->Delete("p");
  }
  TGraphErrors* gMes = new TGraphErrors(xMes.size(),&(xMes[0]),&(yMes[0]),&(xErr[0]),&(yErr[0]));
  gMes->SetLineWidth(1);
  gMes->SetLineColor(kBlue);
  gMes->SetMarkerColor(kBlue);
  gMes->SetMarkerStyle(kFullCircle);
  gMes->SetMarkerSize(0.1);
  gMes->Draw("p");

  ifstream fseeds(Form("seeds/event%09d-Seed.csv",selectedEvent));
  for( std::string line; getline(fseeds, line); ){
    TObjArray* objArray = TString(line).Tokenize(",");
    if (TString(objArray->At(0)->GetName()).Contains("seed_id")) continue;
    double x[3];
    double y[3];
    x[0] = TString(objArray->At( 5)->GetName()).Atof();
    y[0] = TString(objArray->At( 6)->GetName()).Atof();
    x[1] = TString(objArray->At( 8)->GetName()).Atof();
    y[1] = TString(objArray->At( 9)->GetName()).Atof();
    x[2] = TString(objArray->At(11)->GetName()).Atof();
    y[2] = TString(objArray->At(12)->GetName()).Atof();
    TGraph* g = new TGraph(3,x,y);
    g->SetLineWidth(1);
    g->SetLineColor(kRed);
    g->SetMarkerColor(kRed);
    g->SetMarkerStyle(kFullCircle);
    g->SetMarkerSize(0.5);
    g->Draw("lp");
    // TEllipse* bSeed = new TEllipse(bx,by,5);
    // bSeed->SetFillColor(kRed);
    // bSeed->SetLineColor(kRed);
    // bSeed->Draw();
    // TEllipse* mSeed = new TEllipse(mx,my,5);
    // mSeed->SetFillColor(kRed);
    // mSeed->SetLineColor(kRed);
    // mSeed->Draw();
    // TEllipse* tSeed = new TEllipse(tx,ty,5);
    // tSeed->SetFillColor(kRed);
    // tSeed->SetLineColor(kRed);
    // tSeed->Draw();
    objArray->Delete();
  }

}
