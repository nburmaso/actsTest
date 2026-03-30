#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "style.h"
#include "../MyFtdGeo.h"
// for single-track events only
//void analyse_spacepoint_resolution(TString dir = "../build/test02", int selected_layer = 4){
//void analyse_spacepoint_resolution(TString dir = "../build/test02_18", int selected_layer = 4){
//void analyse_spacepoint_resolution(TString dir = "../build/test02_16", int selected_layer = 4){
void analyse_spacepoint_resolution(TString dir = "../build/test", int selected_layer = 1){
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
  vector<vector<int>> vm(nEvents);
  vector<double> vmx;
  vector<double> vmy;  

  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    int layerType = fg.GetLayerType(m_layer_id-2);
    if (layerType == 2) continue;
    if (fabs(m_true_z/10-lz)>3.1) continue;
    vn[m_event_id]++;
    vm[m_event_id].push_back(m_layer_id-2);
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
    if (fabs(tz/10-lz)>0.01) continue;
    vx[event_id] = tx;
    vy[event_id] = ty;
  }

  TFile* fSpacepoints = new TFile(dir + "spacepoints.root");
  TTree* tSpacepoints = (TTree*) fSpacepoints->Get("spacepoints");
  tSpacepoints->Print();
  float sx;
  float sy;
  float sz;
  float st;  
  float varxx;
  float varxy;
  float varyy;
  UInt_t sevent_id;
  ULong64_t sgeometry_id;
  tSpacepoints->SetBranchAddress("x",&sx);
  tSpacepoints->SetBranchAddress("y",&sy);
  tSpacepoints->SetBranchAddress("z",&sz);
  tSpacepoints->SetBranchAddress("t",&varxy);
  tSpacepoints->SetBranchAddress("var_r",&varxx);  
  tSpacepoints->SetBranchAddress("var_z",&varyy);  
  tSpacepoints->SetBranchAddress("geometry_id",&sgeometry_id);
  tSpacepoints->SetBranchAddress("event_id",&sevent_id);
  TH2D* hXY = new TH2D("hXY","",260,-1300,1300,260,-1300,1300);  
  TH2D* hDxDy = new TH2D("hDxDy","",100,-10,10,100,-10,10);  
  TH1D* hDr = new TH1D("hDr","",100,-10,10);  
  TH1D* hDr3 = new TH1D("hDr3","",100,-10,10);  
  TH1D* hDr4 = new TH1D("hDr4","",100,-10,10);  
  TH1D* hDr5 = new TH1D("hDr5","",100,-10,10);
  
  TH1D* hDp = new TH1D("hDp","",1000,-1,1);    
  TH1D* hDp3 = new TH1D("hDp3","",1000,-1,1);    
  TH1D* hDp4 = new TH1D("hDp4","",1000,-1,1);    
  TH1D* hDp5 = new TH1D("hDp5","",1000,-1,1);

  TH1D* hk = new TH1D("hk","",100,-1.e-3,1.e-3);

  TH1D* hdx = new TH1D("hdx","",100,0,1.);
  TH1D* hdy = new TH1D("hdy","",100,0,1.);
  TH1D* hdxdy = new TH1D("hdxdy","",100,0,1.);

  TH1D* hdr = new TH1D("hdr","",100,0,0.8);
  TH1D* hdr3 = new TH1D("hdr3","",100,0,0.8);
  TH1D* hdr4 = new TH1D("hdr4","",100,0,0.8);
  TH1D* hdr5 = new TH1D("hdr5","",100,0,0.8);
  TH1D* hdp = new TH1D("hdp"  ,"",120,0,0.12);
  TH1D* hdp3 = new TH1D("hdp3","",120,0,0.12);
  TH1D* hdp4 = new TH1D("hdp4","",120,0,0.12);
  TH1D* hdp5 = new TH1D("hdp5","",120,0,0.12);

  TH1D* hPull_r = new TH1D("hPull_r","",100,-5,5);
  TH1D* hPull_r3 = new TH1D("hPull_r3","",100,-5,5);
  TH1D* hPull_r4 = new TH1D("hPull_r4","",100,-5,5);
  TH1D* hPull_r5 = new TH1D("hPull_r5","",100,-5,5);
  TH1D* hPull_p = new TH1D("hPull_p"  ,"",120,-5,5);
  TH1D* hPull_p3 = new TH1D("hPull_p3","",120,-5,5);
  TH1D* hPull_p4 = new TH1D("hPull_p4","",120,-5,5);
  TH1D* hPull_p5 = new TH1D("hPull_p5","",120,-5,5);


  for (int is=0;is<tSpacepoints->GetEntries();is++){
    tSpacepoints->GetEntry(is);
    if (fabs(sz/10-lz)>0.1) continue;
    double hx = vx[sevent_id];
    double hy = vy[sevent_id];    
    double dx = sx - hx;
    double dy = sy - hy;
    double sr = sqrt(sx*sx+sy*sy);
    double hr = sqrt(hx*hx+hy*hy);    
    double dr = sr - hr;
    double dp = hr*(atan2(sy,sx)-atan2(hy,hx));

    double ht = atan2(hy,hx);
    double varrr = (hx*hx*varxx + hy*hy*varyy + 2*hx*hy*varxy)/(hr*hr);
    double ha = (ht*hx - hy);
    double hb = (ht*hy + hx);
    double varpp = (hy*hy*varxx + hx*hx*varyy - 2*hx*hy*varxy)/(hr*hr);


    // if (sy<100) continue;
    // if (sx<100) continue;
    hDr->Fill(dr);
    hDp->Fill(dp);
    if (vn[sevent_id]==3) hDr3->Fill(dr);
    if (vn[sevent_id]==4) hDr4->Fill(dr);    
    if (vn[sevent_id]==5) hDr5->Fill(dr);

    if (vn[sevent_id]==3) hDp3->Fill(dp);
    if (vn[sevent_id]==4) hDp4->Fill(dp);    
    if (vn[sevent_id]==5) hDp5->Fill(dp);
    printf("varxx=%e varyy=%e varxy=%e\n",varxx,varyy,varxy);
    hdx->Fill(sqrt(varxx));
    hdy->Fill(sqrt(varyy));    
    hdxdy->Fill(sqrt(fabs(varxy)));    
    hdr->Fill(sqrt(varrr));
    hdp->Fill(sqrt(varpp));    
    if (vn[sevent_id]==3) hdr3->Fill(sqrt(varrr));
    if (vn[sevent_id]==4) hdr4->Fill(sqrt(varrr));
    if (vn[sevent_id]==5) hdr5->Fill(sqrt(varrr));

    if (vn[sevent_id]==3) hdp3->Fill(sqrt(varpp));
    if (vn[sevent_id]==4) hdp4->Fill(sqrt(varpp));
    if (vn[sevent_id]==5) hdp5->Fill(sqrt(varpp));

    double dr_est = sqrt(varrr);
    double dp_est = sqrt(varpp);

    hPull_r->Fill(dr/dr_est);
    hPull_p->Fill(dp/dp_est);
    if (vn[sevent_id]==3) hPull_r3->Fill(dr/dr_est);
    if (vn[sevent_id]==4) hPull_r4->Fill(dr/dr_est);
    if (vn[sevent_id]==5) hPull_r5->Fill(dr/dr_est);

    if (vn[sevent_id]==3) hPull_p3->Fill(dp/dp_est);
    if (vn[sevent_id]==4) hPull_p4->Fill(dp/dp_est);
    if (vn[sevent_id]==5) hPull_p5->Fill(dp/dp_est);

    hXY->Fill(sx,sy);
    hDxDy->Fill(dx,dy);
    if (dp>0.46 && vn[sevent_id]==4) {//printf("event=%d n=%d\n",sevent_id,vn[sevent_id]);
    //if (fabs(dp)<0.1 && vn[sevent_id]==4) {
    //if (fabs(dp)>0.) {
      printf("event=%d n=%d %f ",sevent_id,vn[sevent_id],sqrt(st));
      for (auto k : vm[sevent_id]){
        printf("%d ",k);
      }
      printf("\n");
    }
//    if (fabs(dp)>0.4) printf("event=%d n=%d\n",sevent_id,vn[sevent_id]);
  }

  TH1D* hN = new TH1D("hN","",10,0,10);
  for (int ev=0;ev<vn.size();ev++){
    hN->Fill(vn[ev]);
  }

  new TCanvas;
  hdx->Draw();
  
  new TCanvas;
  hdy->Draw();
  
  new TCanvas;
  hdxdy->Draw();
  
  new TCanvas;
  hdr->Draw();
  hdr3->SetLineColor(kRed);
  hdr4->SetLineColor(kMagenta);
  hdr5->SetLineColor(kBlue);
  hdr3->Draw("same");
  hdr4->Draw("same");
  hdr5->Draw("same");
  
  new TCanvas;
  hdp->Draw();
  hdp3->SetLineColor(kRed);
  hdp4->SetLineColor(kMagenta);
  hdp5->SetLineColor(kBlue);
  hdp3->Draw("same");  
  hdp4->Draw("same");  
  hdp5->Draw("same");

  new TCanvas;
  hPull_r->Draw();
  hPull_r3->SetLineColor(kRed);
  hPull_r4->SetLineColor(kMagenta);
  hPull_r5->SetLineColor(kBlue);
  hPull_r3->Draw("same");
  hPull_r4->Draw("same");
  hPull_r5->Draw("same");
  hPull_r5->Fit("gaus");
  
  new TCanvas;
  hPull_p->Draw();
  hPull_p3->SetLineColor(kRed);
  hPull_p4->SetLineColor(kMagenta);
  hPull_p5->SetLineColor(kBlue);
  hPull_p3->Draw("same");
  hPull_p4->Draw("same");
  hPull_p5->Draw("same");
  hPull_p5->Fit("gaus");
}
