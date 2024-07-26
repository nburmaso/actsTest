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
#include "ActsFatras/EventData/Barcode.hpp"

using namespace std;

//void draw(TString dir = "mpd"){
void draw(TString dir = "mpd_notpc"){
//void draw(TString dir = "notpc_pi_16"){
  dir.Append("/");
  
  int selectedEvent = 55;
  UInt_t part_event_id;
  auto part_particle_id   = new vector<unsigned long>; 
  auto particle_type = new vector<int>;
  auto vz = new vector<float>;
  auto px = new vector<float>;
  auto py = new vector<float>;
  auto pz = new vector<float>;
  auto pt = new vector<float>;

  TFile* fPart = new TFile(dir + "particles.root");
  TTree* tPart = (TTree*) fPart->Get("particles");
  tPart->SetBranchAddress("event_id",&part_event_id);
  tPart->SetBranchAddress("particle_id",&part_particle_id);
  tPart->SetBranchAddress("particle_type",&particle_type);
  tPart->SetBranchAddress("vz",&vz);
  tPart->SetBranchAddress("px",&px);
  tPart->SetBranchAddress("py",&py);
  tPart->SetBranchAddress("pz",&pz);  
  tPart->SetBranchAddress("pt",&pt);  
  tPart->GetEntry(selectedEvent);

  // for (int i=0;i<100;i++){
  //   tPart->GetEntry(i);
  //   printf("%d %f\n", i, pt->at(98));
  // }
  //return;

  TFile* fSeeds = new TFile(dir + "seeds.root");
  TTree* tSeeds = (TTree*) fSeeds->Get("seeds");
  UInt_t seed_event_id;
  ULong64_t measurement_id_1;
  ULong64_t measurement_id_2;
  ULong64_t measurement_id_3;
  tSeeds->SetBranchAddress("event_id", &seed_event_id);
  tSeeds->SetBranchAddress("measurement_id_1", &measurement_id_1);
  tSeeds->SetBranchAddress("measurement_id_2", &measurement_id_2);
  tSeeds->SetBranchAddress("measurement_id_3", &measurement_id_3);

  TFile* fMes = new TFile(dir+ "measurements.root");
  TTree* tMes = (TTree*) fMes->Get("vol1");
  // tMes->Print();
  float true_x;
  float true_y;
  float true_z;
  int event_nr;
  tMes->SetBranchAddress("true_x",&true_x);
  tMes->SetBranchAddress("true_y",&true_y);
  tMes->SetBranchAddress("true_z",&true_z);
  tMes->SetBranchAddress("event_nr",&event_nr);

  std::map<int,int> mapFirstMeaurementIndex;
  int previous_event_nr = -1;
  for (int iMes=0;iMes<tMes->GetEntries();iMes++){
    tMes->GetEntry(iMes);
    if (event_nr==previous_event_nr) continue;
    previous_event_nr = event_nr;
    mapFirstMeaurementIndex[event_nr] = iMes;
  }

  for (int iSeeds=0;iSeeds<tSeeds->GetEntries();iSeeds++){
    tSeeds->GetEntry(iSeeds);
    if (seed_event_id!=selectedEvent) continue;

    printf("%d %llu %llu %llu\n", seed_event_id, measurement_id_1, measurement_id_2, measurement_id_3);
    int shift = mapFirstMeaurementIndex[seed_event_id];
    tMes->GetEntry(shift+measurement_id_1);
    printf("%f\n", true_z);
    tMes->GetEntry(shift+measurement_id_2);
    printf("%f\n", true_z);
    tMes->GetEntry(shift+measurement_id_3);
    printf("%f\n", true_z);
  }


  TFile* fHits = new TFile(dir + "hits.root");
  TTree* tHits = (TTree*) fHits->Get("hits");
  float tz = 0;
  float tx = 0;
  float ty = 0;
  UInt_t event_id = 0;
  ULong64_t particle_id = 0;
  ULong64_t geometry_id = 0;
  tHits->SetBranchAddress("event_id",&event_id);
  tHits->SetBranchAddress("particle_id",&particle_id);
  tHits->SetBranchAddress("tx",&tx);
  tHits->SetBranchAddress("ty",&ty);
  tHits->SetBranchAddress("tz",&tz);
  tHits->BuildIndex("event_id","int(tx*1048576)");

  for (int iMes=0;iMes<tMes->GetEntries();iMes++){
    tMes->GetEntry(iMes);
    if (event_nr!=selectedEvent) continue;
    // no hit index in measurement tree. using workaround
    tHits->GetEntryWithIndex(event_nr, int(true_x*1048576));
    if (abs(true_y-ty)>1e-10) {
      printf("WARNING Wrong hit index\n");
    }
  }

  TH1D* hZ = new TH1D("hZ","",100,2000,3100);
  const int nMax = 1000;
  double x[nMax];
  double y[nMax];
  int n=0;
  printf("%lld\n",tHits->GetEntries());

  TH1D* hPartZ = new TH1D("PartZ","",100,0,3100);

  std::map<ULong64_t,TGraph*> mapGraphHitsXY;
  for (int entry=0;entry< tHits->GetEntries();entry++){
    tHits->GetEntry(entry);
    //printf("%d",event_id);
    if (event_id!=selectedEvent) continue;
    hZ->Fill(tz);
    printf("%llu\n",particle_id);
    if (auto graph = mapGraphHitsXY.find(particle_id); graph == mapGraphHitsXY.end()) mapGraphHitsXY[particle_id] = new TGraph();
    mapGraphHitsXY[particle_id]->AddPoint(tx,ty);

    ActsFatras::Barcode barcode(particle_id);
    int partIndex = barcode.particle()-1;
    float part_vz = vz->at(partIndex);
    float part_pt = pt->at(partIndex);
    printf("%f\n",part_vz);
    hPartZ->Fill(part_vz);
    if (part_pt>0.1) {
      mapGraphHitsXY[particle_id]->SetLineColor(kGreen);
      mapGraphHitsXY[particle_id]->SetMarkerColor(kGreen);
    }
    if (abs(part_vz)<1000) continue;
    // secondary
    mapGraphHitsXY[particle_id]->SetLineColor(kGray);
    mapGraphHitsXY[particle_id]->SetMarkerColor(kGray);
  }

  new TCanvas;
  hPartZ->Draw();

  new TCanvas;
  hZ->Draw();
  new TCanvas("chits","chits",900,920);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.01);

  TH1F* hFrame = gPad->DrawFrame(-1300,-1300,1300,1300);
  hFrame->SetTitle("; mm;");
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

  for (int iSeeds=0;iSeeds<tSeeds->GetEntries();iSeeds++){
    tSeeds->GetEntry(iSeeds);
    if (seed_event_id!=selectedEvent) continue;
    double x[3];
    double y[3];
    int shift = mapFirstMeaurementIndex[seed_event_id];
    tMes->GetEntry(shift+measurement_id_1);
    x[0] = true_x;
    y[0] = true_y;
    tMes->GetEntry(shift+measurement_id_2);
    x[1] = true_x;
    y[1] = true_y; 
    tMes->GetEntry(shift+measurement_id_3);
    x[2] = true_x; 
    y[2] = true_y; 
    TGraph* g = new TGraph(3,x,y);
    g->SetLineWidth(1);
    g->SetLineColor(kRed);
    g->SetMarkerColor(kRed);
    g->SetMarkerStyle(kFullCircle);
    g->SetMarkerSize(0.5);
    g->Draw("lp");
  }
  gPad->Print(dir + "display.png");
}
