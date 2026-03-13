#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "style.h"
#include "../MyFtdGeo.h"
// for single-track events only
//void analyse_spacepoint_resolution(TString dir = "../build/test02", int selected_layer = 4){
void analyse_spacepoint_resolution(TString dir = "../build/test02", int selected_layer = 4){
  dir.Append("/");
  MyFtdGeo fg;
  double lz = fg.GetLayerPositions()[selected_layer];

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
  
  vector<int> vn(nEvents,0);
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    int layerType = fg.GetLayerType(m_layer_id-2);
    if (layerType == 2) continue;
    if (fabs(m_true_z/10-lz)>3.1) continue;
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
    if (fabs(tz/10-lz)>1) continue;
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
  TH1D* hDr3 = new TH1D("hDr3","",100,-10,10);  
  TH1D* hDr4 = new TH1D("hDr4","",100,-10,10);  
  TH1D* hDr5 = new TH1D("hDr5","",100,-10,10);
  
  TH1D* hDp = new TH1D("hDp","",100,-2,2);    
  TH1D* hDp3 = new TH1D("hDp3","",100,-2,2);    
  TH1D* hDp4 = new TH1D("hDp4","",100,-2,2);    
  TH1D* hDp5 = new TH1D("hDp5","",100,-2,2);    
  for (int is=0;is<tSpacepoints->GetEntries();is++){
    tSpacepoints->GetEntry(is);
    if (fabs(sz/10-lz)>1) continue;
    double hx = vx[sevent_id];
    double hy = vy[sevent_id];    
    double dx = sx - hx;
    double dy = sy - hy;
    double sr = sqrt(sx*sx+sy*sy);
    double hr = sqrt(hx*hx+hy*hy);    
    double dr = sr - hr;
    double dp = hr*(atan2(sy,sx)-atan2(hy,hx));
    // if (fabs(sy)>100) continue;
    // if (sx<0) continue;
    hDr->Fill(dr);
    hDp->Fill(dp);
    if (vn[sevent_id]==3) hDr3->Fill(dr);
    if (vn[sevent_id]==4) hDr4->Fill(dr);    
    if (vn[sevent_id]==5) hDr5->Fill(dr);

    if (vn[sevent_id]==3) hDp3->Fill(dp);
    if (vn[sevent_id]==4) hDp4->Fill(dp);    
    if (vn[sevent_id]==5) hDp5->Fill(dp);

    hXY->Fill(sx,sy);
    hDxDy->Fill(dx,dy);
    if (fabs(dp)<0.1) printf("event=%d n=%d\n",sevent_id,vn[sevent_id]);
  }

  TH1D* hN = new TH1D("hN","",10,0,10);
  for (int ev=0;ev<vn.size();ev++){
    hN->Fill(vn[ev]);
  }
  new TCanvas;
  hN->Draw();

  hXY->Draw();
  new TCanvas;
  hDxDy->Draw();
  new TCanvas;
  hDr->Draw();
  hDr3->SetLineColor(kRed);
  hDr4->SetLineColor(kMagenta);
  hDr5->SetLineColor(kBlue);  
  hDr3->Draw("same");
  hDr4->Draw("same");
  hDr5->Draw("same");
  new TCanvas;
  hDp->Draw();
  hDp3->SetLineColor(kRed);
  hDp4->SetLineColor(kMagenta);
  hDp5->SetLineColor(kBlue);  
  hDp3->Draw("same");
  hDp4->Draw("same");
  hDp5->Draw("same");
  hDp4->Fit("gaus","","",-0.2,0.2);
}
