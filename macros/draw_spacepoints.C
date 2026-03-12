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
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "TBox.h"
#include "TPolyLine.h"
#include "Math/Point2D.h"
#include "Math/Vector2D.h"
using namespace std;

void draw_spacepoints(int selected_event = 0, int selected_layer = 0, int zoom = 0){
  TString dir = "../build/test";
  dir.Append("/");

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

  float sx;
  float sy;
  float sz;
  UInt_t sevent_id;
  ULong64_t sgeometry_id;
  TFile* fSpacepoints = new TFile(dir + "spacepoints.root");
  TTree* tSpacepoints = (TTree*) fSpacepoints->Get("spacepoints");
  tSpacepoints->Print();
  // return;

  tSpacepoints->SetBranchAddress("x",&sx);
  tSpacepoints->SetBranchAddress("y",&sy);
  tSpacepoints->SetBranchAddress("z",&sz);
  tSpacepoints->SetBranchAddress("geometry_id",&sgeometry_id);
  tSpacepoints->SetBranchAddress("event_id",&sevent_id);

  // new TCanvas("cmeas","cmeas",1320,1320);
  new TCanvas("cmeas","cmeas",850,850);
  gPad->SetRightMargin(0.005);
  gPad->SetTopMargin(0.005);
  gPad->SetLeftMargin(0.1);
  gPad->SetBottomMargin(0.07);

  TH1F* hFrame = gPad->DrawFrame(-150, -150, 150, 150);
  // if (zoom==1) hFrame = gPad->DrawFrame(50, -25, 95, 20);

  hFrame->GetYaxis()->SetTitleOffset(1.5);
  hFrame->SetTitle(";x (cm); y (cm)"); 

  TGraph* gHits = new TGraph();
  for (int ih=0;ih<tHits->GetEntries();ih++){
    tHits->GetEntry(ih);
    if (selected_event>=0 && event_id!=selected_event) continue;
    // if (sz<2069 || sz>2071) continue;
    // if (tz<2339 || tz>2341) continue;
    gHits->AddPoint(tx/10.,ty/10.);
  }

  TGraph* gSP = new TGraph();
  for (int is=0;is<tSpacepoints->GetEntries();is++){
      tSpacepoints->GetEntry(is);
      if (selected_event>=0 && sevent_id!=selected_event) continue;
      auto geoId = Acts::GeometryIdentifier(sgeometry_id);
      int layerActs = geoId.layer();
      int surfaceActs = geoId.sensitive();
      if (surfaceActs==1) continue;
      //printf("%d\n",layerActs);
//      if (sz<2069 || sz>2071) continue;
//      if (sz<2339 || sz>2341) continue;
//      if (sz<2779 || sz>2781) continue;

      cout << geoId << endl;
      // int layerId = (layerActs - 2); // Mpd notation
      // if (layerId < selected_layer || layerId >= selected_layer + 7) continue;
      // printf("%d %f %f %f\n", layerId, sx, sy, sz);
      gSP->AddPoint(sx/10.,sy/10.);
    }
  gHits->SetMarkerColor(kBlue);
  gHits->SetMarkerStyle(kOpenCircle);
  gHits->SetMarkerSize(1.0);
  gHits->Draw("p");
  gSP->SetMarkerColor(kMagenta);
  gSP->SetMarkerStyle(kFullCircle);
  gSP->SetMarkerSize(1.0);
  gSP->Draw("p");
  // gPad->Print(Form("event_%d_station_%d.png",selected_event,selected_layer));
}
