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
#include "../MyFtdGeo.h"
using namespace std;

const int shift = 2; //  isroc = 1;  isframe = 1; without fake pre-layer

// set layer corresponding to fake surface in the middle of the station
void draw_straw_station(int selected_event = 0, int selected_layer = 4, bool draw_spacepoints = 1, bool zoom = 0){
  TString dir = "../build/test";
  dir.Append("/");

  MyFtdGeo fg;
  double incl = fg.GetTubeIncl();
  double rmin = fg.GetLayerRMin(selected_layer);
  double rmax = fg.GetLayerRMax(selected_layer);
  double dr = (rmax-rmin)/2.;
  double rc = (rmin+rmax)/2.;
  double lz = fg.GetLayerPositions()[selected_layer];
  const int nLayersPerStation = fg.GetNLayersPerStation();

  // new TCanvas("cmeas","cmeas",1320,1320);
  new TCanvas("cmeas","cmeas",850,850);
  gPad->SetRightMargin(0.005);
  gPad->SetTopMargin(0.005);
  gPad->SetLeftMargin(0.1);
  gPad->SetBottomMargin(0.07);
  TH1F* hFrame = gPad->DrawFrame(-rmax, -rmax, rmax, rmax);
  hFrame->GetYaxis()->SetTitleOffset(1.5);
  hFrame->SetTitle(";x (cm); y (cm)"); 

  TBox* b = new TBox();
  b->SetLineColor(kBlue-10);
  b->SetLineWidth(1);
  b->SetLineStyle(1);
  b->SetFillColor(kRed);
  b->SetFillStyle(0);
  b->Draw();

  using std::numbers::pi;

  std::map<int, TPolyLine*> mapPolyLines[nLayersPerStation];
  for (int layer=selected_layer-nLayersPerStation/2; layer<=selected_layer+nLayersPerStation/2; layer++){
    int layerType = fg.GetLayerType(layer);
    if (layerType == 2) continue;
    bool back = layer>selected_layer;
    auto color = back ? kOrange : kOrange+7;
    if (layerType == 5) color = back ? kGreen-7 : kGreen+1;
    if (layerType == 6) color = back ? kAzure+6 : kAzure;
    int nStraws = fg.GetLayerNumberOfTubes(layer);
    double strawHalfLength = dr/cos(fg.GetLayerStereoAngle(layer));
    for (int strawId=0; strawId < nStraws; strawId++){
      double xc = fg.GetTubeCenterX(layer, strawId);
      double yc = fg.GetTubeCenterY(layer, strawId);
      ROOT::Math::XYPoint pc(xc, yc);
      ROOT::Math::Polar2DVector vc1(strawHalfLength, fg.GetTubeRotationAngle(layer, strawId));
      ROOT::Math::XYPoint prmax = pc + vc1;
      ROOT::Math::XYPoint prmin = pc - vc1;
      double x[] = {prmin.x(), prmax.x()};
      double y[] = {prmin.y(), prmax.y()};
      TPolyLine *pline = new TPolyLine(2, x, y);
      pline->SetFillColor(color);
      pline->SetLineColor(color);
      pline->SetLineWidth(1);
      // pline->Draw("f");
      // pline->Draw();
      mapPolyLines[layer%nLayersPerStation].emplace(strawId, pline);
    }
  }

  TFile* fMeas = new TFile(dir+"measurements.root");
  TTree* tMeas = (TTree*) fMeas->Get("measurements");
  // tMeas->Print();
  int32_t m_event_id;
  int32_t m_volume_id;
  int32_t m_layer_id;
  int32_t m_surface_id;
  float m_true_loc0;
  float m_true_z;
  tMeas->SetBranchAddress("event_nr",&m_event_id);
  tMeas->SetBranchAddress("volume_id",&m_volume_id);
  tMeas->SetBranchAddress("layer_id",&m_layer_id);
  tMeas->SetBranchAddress("surface_id",&m_surface_id);
  tMeas->SetBranchAddress("true_loc0",&m_true_loc0);
  tMeas->SetBranchAddress("true_z",&m_true_z);
  
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if (selected_event>=0 && m_event_id!=selected_event) continue;
    int layer = m_layer_id-shift;
    int layerType = fg.GetLayerType(layer);
    if (layerType == 2) continue;
    if (fabs(m_true_z/10-lz)>3) continue;
    printf("%d %d\n",m_layer_id,m_surface_id);
    auto pline = mapPolyLines[layer%nLayersPerStation][m_surface_id-1];
    pline->SetLineWidth(2);
    pline->Draw();
  }


  TEllipse* elDot = new TEllipse(0,0, 2);
  TEllipse* elMin = new TEllipse(0,0,rmin);
  TEllipse* elMax = new TEllipse(0,0,rmax);
  elDot->SetFillStyle(0);
  elMin->SetFillStyle(0);
  elMax->SetFillStyle(0);
  elMax->Draw();
  elMin->Draw();
  elDot->Draw();

  TFile* fHits = new TFile(dir + "hits.root");
  TTree* tHits = (TTree*) fHits->Get("hits");
  float tz = 0;
  float tx = 0;
  float ty = 0;
  UInt_t event_id = 0;
  ULong64_t geometry_id = 0;
  tHits->SetBranchAddress("event_id",&event_id);
  tHits->SetBranchAddress("geometry_id",&geometry_id);
  tHits->SetBranchAddress("tx",&tx);
  tHits->SetBranchAddress("ty",&ty);
  tHits->SetBranchAddress("tz",&tz);
  
  TGraph* g = new TGraph();
  TGraph* gh = new TGraph();  
  for (int ih=0;ih<tHits->GetEntries();ih++){
    tHits->GetEntry(ih);
    if (selected_event>=0 && event_id!=selected_event) continue;
    auto geoId = Acts::GeometryIdentifier(geometry_id);
    int layer = geoId.layer()-shift;      
    if (fabs(layer-selected_layer)>nLayersPerStation/2) continue;
    g->AddPoint(tx/10,ty/10);
    if (fabs(tz/10-lz)<0.01) gh->AddPoint(tx/10,ty/10);
  }

//  g->SetMarkerColor(color);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerSize(1.);
  g->Draw("p");
  gh->SetMarkerStyle(kFullCircle);
  gh->SetMarkerSize(1.);
  gh->SetMarkerColor(kGreen);
  gh->Draw("p");

  float sx;
  float sy;
  float sz;
  UInt_t sevent_id;
  ULong64_t sgeometry_id;

  TFile* fSpacepoints = new TFile(dir + "spacepoints.root");
  TTree* tSpacepoints = (TTree*) fSpacepoints->Get("spacepoints");

  tSpacepoints->SetBranchAddress("x",&sx);
  tSpacepoints->SetBranchAddress("y",&sy);
  tSpacepoints->SetBranchAddress("z",&sz);
  tSpacepoints->SetBranchAddress("geometry_id",&sgeometry_id);
  tSpacepoints->SetBranchAddress("event_id",&sevent_id);

  TGraph* gSP = new TGraph();
  for (int is=0;is<tSpacepoints->GetEntries();is++){
      tSpacepoints->GetEntry(is);
      if (selected_event>=0 && sevent_id!=selected_event) continue;
      auto geoId = Acts::GeometryIdentifier(sgeometry_id);
      int layer = geoId.layer()-shift;      
      if (fabs(layer-selected_layer)>nLayersPerStation/2) continue;
      gSP->AddPoint(sx/10.,sy/10.);
  }
  gSP->SetMarkerColor(kMagenta);
  gSP->SetMarkerStyle(kFullCircle);
  gSP->SetMarkerSize(0.7);
  gSP->Draw("p");

}
