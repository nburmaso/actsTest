#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "style.h"
#include "../MyFtdGeo.h"
// for single-track events only
void analyse_spacepoint_resolution(TString dir = "../build/test02", int selected_layer = 0){
  dir.Append("/");
  MyFtdGeo fg;
  double lz = fg.GetLayerPositions()[selected_layer+2];

  TFile* fPart = new TFile(dir + "particles.root");
  TTree* tPart = (TTree*) fPart->Get("particles");
  int nEvents = tPart->GetEntries();

  TFile* fMeas = new TFile(dir + "measurements.root");
  TTree* tMeas = (TTree*) fMeas->Get("measurements");
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
  
  vector<double> vn(nEvents,0);
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if (fabs(m_true_z/10-lz)>2.1) continue;
    vn[m_event_id]++;
  }

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

  vector<double> vx(nEvents,0);
  vector<double> vy(nEvents,0);

  for (int ih=0;ih<tHits->GetEntries();ih++){
    tHits->GetEntry(ih);
    if (fabs(tz/10-lz)>2) continue;
    vx[event_id] = tx;
    vy[event_id] = ty;
  }

  TFile* fSpacepoints = new TFile(dir + "spacepoints.root");
  TTree* tSpacepoints = (TTree*) fSpacepoints->Get("spacepoints");
  float sx;
  float sy;
  float sz;
  UInt_t sevent_id;
  ULong64_t sgeometry_id;
  tSpacepoints->SetBranchAddress("x",&sx);
  tSpacepoints->SetBranchAddress("y",&sy);
  tSpacepoints->SetBranchAddress("z",&sz);
  tSpacepoints->SetBranchAddress("geometry_id",&sgeometry_id);
  tSpacepoints->SetBranchAddress("event_id",&sevent_id);

  TH2D* hXY = new TH2D("hXY","",260,-1300,1300,260,-1300,1300);  
  TH2D* hDxDy = new TH2D("hDxDy","",100,-10,10,100,-10,10);  
  TH1D* hDr = new TH1D("hDr","",100,-10,10);  
  for (int is=0;is<tSpacepoints->GetEntries();is++){
    tSpacepoints->GetEntry(is);
    if (fabs(sz/10-lz)>2) continue;
    double hx = vx[sevent_id];
    double hy = vy[sevent_id];    
    double dx = sx - hx;
    double dy = sy - hy;
    double sr = sqrt(sx*sx+sy*sy);
    double hr = sqrt(hx*hx+hy*hy);    
    double dr = sr - hr;
    hDr->Fill(dr);
    if (fabs(sy)>100) continue;
    if (sx<0) continue;
    hXY->Fill(sx,sy);
    hDxDy->Fill(dx,dy);
  }

  hXY->Draw();
  new TCanvas;
  hDxDy->Draw();
  new TCanvas;
  hDr->Draw();

  new TCanvas;
  TH1D* hN = new TH1D("hN","",10,0,10);
  for (int ev=0;ev<vn.size();ev++){
    hN->Fill(vn[ev]);
  }
  hN->Draw();
}
