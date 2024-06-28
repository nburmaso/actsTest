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

#include "style.h"

using namespace std;

// data from RootTrackSummaryReader.hpp
  /// the event number
  uint32_t m_eventNr{0};

  /// The entry numbers for accessing events in increased order (there could be
  /// multiple entries corresponding to one event number)
  std::vector<long long> m_entryNumbers = {};

  /// The number of states
  std::vector<unsigned int>* m_nStates = new std::vector<unsigned int>;
  /// The number of measurements
  std::vector<unsigned int>* m_nMeasurements = new std::vector<unsigned int>;
  /// The number of outliers
  std::vector<unsigned int>* m_nOutliers = new std::vector<unsigned int>;
  /// The number of holes
  std::vector<unsigned int>* m_nHoles = new std::vector<unsigned int>;
  /// The total chi2
  std::vector<float>* m_chi2Sum = new std::vector<float>;
  /// The number of ndf of the measurements+outliers
  std::vector<unsigned int>* m_NDF = new std::vector<unsigned int>;
  /// The chi2 on all measurement states
  std::vector<std::vector<double>>* m_measurementChi2 =
      new std::vector<std::vector<double>>;
  /// The chi2 on all outlier states
  std::vector<std::vector<double>>* m_outlierChi2 =
      new std::vector<std::vector<double>>;
  /// The volume id of the measurements
  std::vector<std::vector<std::uint32_t>>* m_measurementVolume =
      new std::vector<std::vector<std::uint32_t>>;
  /// The layer id of the measurements
  std::vector<std::vector<std::uint32_t>>* m_measurementLayer =
      new std::vector<std::vector<std::uint32_t>>;
  /// The volume id of the outliers
  std::vector<std::vector<std::uint32_t>>* m_outlierVolume =
      new std::vector<std::vector<std::uint32_t>>;
  /// The layer id of the outliers
  std::vector<std::vector<std::uint32_t>>* m_outlierLayer =
      new std::vector<std::vector<std::uint32_t>>;

  // The majority truth particle info
  /// The number of hits from majority particle
  std::vector<unsigned int>* m_nMajorityHits = new std::vector<unsigned int>;
  /// The particle Id of the majority particle
  std::vector<uint64_t>* m_majorityParticleId = new std::vector<uint64_t>;
  /// Charge of majority particle
  std::vector<int>* m_t_charge = new std::vector<int>;
  /// Time of majority particle
  std::vector<float>* m_t_time = new std::vector<float>;
  /// Vertex x positions of majority particle
  std::vector<float>* m_t_vx = new std::vector<float>;
  /// Vertex y positions of majority particle
  std::vector<float>* m_t_vy = new std::vector<float>;
  /// Vertex z positions of majority particle
  std::vector<float>* m_t_vz = new std::vector<float>;
  /// Initial momenta px of majority particle
  std::vector<float>* m_t_px = new std::vector<float>;
  /// Initial momenta py of majority particle
  std::vector<float>* m_t_py = new std::vector<float>;
  /// Initial momenta pz of majority particle
  std::vector<float>* m_t_pz = new std::vector<float>;
  /// Initial momenta theta of majority particle
  std::vector<float>* m_t_theta = new std::vector<float>;
  /// Initial momenta phi of majority particle
  std::vector<float>* m_t_phi = new std::vector<float>;
  /// Initial momenta pT of majority particle
  std::vector<float>* m_t_pT = new std::vector<float>;
  /// Initial momenta eta of majority particle
  std::vector<float>* m_t_eta = new std::vector<float>;

  /// If the track has fitted parameter
  std::vector<bool>* m_hasFittedParams = new std::vector<bool>;
  /// Fitted parameters eBoundLoc0 of track
  std::vector<float>* m_eLOC0_fit = new std::vector<float>;
  /// Fitted parameters eBoundLoc1 of track
  std::vector<float>* m_eLOC1_fit = new std::vector<float>;
  /// Fitted parameters ePHI of track
  std::vector<float>* m_ePHI_fit = new std::vector<float>;
  /// Fitted parameters eTHETA of track
  std::vector<float>* m_eTHETA_fit = new std::vector<float>;
  /// Fitted parameters eQOP of track
  std::vector<float>* m_eQOP_fit = new std::vector<float>;
  /// Fitted parameters eT of track
  std::vector<float>* m_eT_fit = new std::vector<float>;
  /// Fitted parameters eLOC err of track
  std::vector<float>* m_err_eLOC0_fit = new std::vector<float>;
  /// Fitted parameters eBoundLoc1 err of track
  std::vector<float>* m_err_eLOC1_fit = new std::vector<float>;
  /// Fitted parameters ePHI err of track
  std::vector<float>* m_err_ePHI_fit = new std::vector<float>;
  /// Fitted parameters eTHETA err of track
  std::vector<float>* m_err_eTHETA_fit = new std::vector<float>;
  /// Fitted parameters eQOP err of track
  std::vector<float>* m_err_eQOP_fit = new std::vector<float>;
  /// Fitted parameters eT err of track
  std::vector<float>* m_err_eT_fit = new std::vector<float>;

  /// Manually added pulls
  std::vector<float>* m_pull_eLOC0_fit = new std::vector<float>;
  std::vector<float>* m_pull_eLOC1_fit = new std::vector<float>;
  std::vector<float>* m_pull_ePHI_fit = new std::vector<float>;
  std::vector<float>* m_pull_eTHETA_fit = new std::vector<float>;
  std::vector<float>* m_pull_eQOP_fit = new std::vector<float>;
  std::vector<float>* m_pull_eT_fit = new std::vector<float>;

void SetBranchAddresses(TTree* m_inputChain){
  // Set the branches
  m_inputChain->SetBranchAddress("event_nr", &m_eventNr);

  // These info is not really stored in the event store, but still read in
  m_inputChain->SetBranchAddress("nStates", &m_nStates);
  m_inputChain->SetBranchAddress("nMeasurements", &m_nMeasurements);
  m_inputChain->SetBranchAddress("nOutliers", &m_nOutliers);
  m_inputChain->SetBranchAddress("nHoles", &m_nHoles);
  m_inputChain->SetBranchAddress("chi2Sum", &m_chi2Sum);
  m_inputChain->SetBranchAddress("NDF", &m_NDF);
  m_inputChain->SetBranchAddress("measurementChi2", &m_measurementChi2);
  m_inputChain->SetBranchAddress("outlierChi2", &m_outlierChi2);
  m_inputChain->SetBranchAddress("measurementVolume", &m_measurementVolume);
  m_inputChain->SetBranchAddress("measurementLayer", &m_measurementLayer);
  m_inputChain->SetBranchAddress("outlierVolume", &m_outlierVolume);
  m_inputChain->SetBranchAddress("outlierLayer", &m_outlierLayer);

  m_inputChain->SetBranchAddress("majorityParticleId", &m_majorityParticleId);
  m_inputChain->SetBranchAddress("nMajorityHits", &m_nMajorityHits);
  m_inputChain->SetBranchAddress("t_charge", &m_t_charge);
  m_inputChain->SetBranchAddress("t_time", &m_t_time);
  m_inputChain->SetBranchAddress("t_vx", &m_t_vx);
  m_inputChain->SetBranchAddress("t_vy", &m_t_vy);
  m_inputChain->SetBranchAddress("t_vz", &m_t_vz);
  m_inputChain->SetBranchAddress("t_px", &m_t_px);
  m_inputChain->SetBranchAddress("t_py", &m_t_py);
  m_inputChain->SetBranchAddress("t_pz", &m_t_pz);
  m_inputChain->SetBranchAddress("t_theta", &m_t_theta);
  m_inputChain->SetBranchAddress("t_phi", &m_t_phi);
  m_inputChain->SetBranchAddress("t_eta", &m_t_eta);
  m_inputChain->SetBranchAddress("t_pT", &m_t_pT);

  m_inputChain->SetBranchAddress("hasFittedParams", &m_hasFittedParams);
  m_inputChain->SetBranchAddress("eLOC0_fit", &m_eLOC0_fit);
  m_inputChain->SetBranchAddress("eLOC1_fit", &m_eLOC1_fit);
  m_inputChain->SetBranchAddress("ePHI_fit", &m_ePHI_fit);
  m_inputChain->SetBranchAddress("eTHETA_fit", &m_eTHETA_fit);
  m_inputChain->SetBranchAddress("eQOP_fit", &m_eQOP_fit);
  m_inputChain->SetBranchAddress("eT_fit", &m_eT_fit);
  m_inputChain->SetBranchAddress("err_eLOC0_fit", &m_err_eLOC0_fit);
  m_inputChain->SetBranchAddress("err_eLOC1_fit", &m_err_eLOC1_fit);
  m_inputChain->SetBranchAddress("err_ePHI_fit", &m_err_ePHI_fit);
  m_inputChain->SetBranchAddress("err_eTHETA_fit", &m_err_eTHETA_fit);
  m_inputChain->SetBranchAddress("err_eQOP_fit", &m_err_eQOP_fit);
  m_inputChain->SetBranchAddress("err_eT_fit", &m_err_eT_fit);

  // added manually
  m_inputChain->SetBranchAddress("pull_eLOC0_fit", &m_pull_eLOC0_fit);
  m_inputChain->SetBranchAddress("pull_eLOC1_fit", &m_pull_eLOC1_fit);
  m_inputChain->SetBranchAddress("pull_ePHI_fit", &m_pull_ePHI_fit);
  m_inputChain->SetBranchAddress("pull_eTHETA_fit", &m_pull_eTHETA_fit);
  m_inputChain->SetBranchAddress("pull_eQOP_fit", &m_pull_eQOP_fit);
  m_inputChain->SetBranchAddress("pull_eT_fit", &m_pull_eT_fit);
}

void analyse_summary(){
  TFile* f = new TFile("build/tracksummary.root");
  TTree* t = (TTree*) f->Get("tracksummary");
  SetBranchAddresses(t);
  int nEvents = t->GetEntries();

  TH1D* h_eQOP_fit = new TH1D("h_eQOP_fit","",200,0,2);
  TH1D* h_eLOC0_fit = new TH1D("h_eLOC0_fit","",200,-90,90);
  TH1D* h_eLOC1_fit = new TH1D("h_eLOC1_fit","",200,-90,90);
  TH1D* h_pull_eQOP_fit = new TH1D("h_pull_eQOP_fit","",200,-10,10);
  TH1D* h_pull_eLOC0_fit = new TH1D("h_pull_eLOC0_fit","",200,-10,10);
  TH1D* h_pull_eLOC1_fit = new TH1D("h_pull_eLOC1_fit","",200,-10,10);

  for (int ev=0; ev<nEvents; ev++){
    t->GetEntry(ev);
    for (int tr=0; tr<m_majorityParticleId->size(); tr++){
      printf("%f\n",m_eQOP_fit->at(tr));
      h_eQOP_fit->Fill(m_eQOP_fit->at(tr));
      h_eLOC0_fit->Fill(m_eLOC0_fit->at(tr));
      h_eLOC1_fit->Fill(m_eLOC1_fit->at(tr));
      h_pull_eQOP_fit->Fill(m_pull_eQOP_fit->at(tr));
      h_pull_eLOC0_fit->Fill(m_pull_eLOC0_fit->at(tr));
      h_pull_eLOC1_fit->Fill(m_pull_eLOC1_fit->at(tr));
    }
  }

  new TCanvas;
  h_eLOC0_fit->Draw();

  new TCanvas;
  h_eLOC1_fit->Draw();

  new TCanvas;
  h_pull_eQOP_fit->Draw();

  new TCanvas;
  h_pull_eLOC0_fit->Draw();

  new TCanvas;
  h_pull_eLOC1_fit->Draw();

  gStyle->SetOptFit(1);

  new TCanvas;
  SetHisto(h_eQOP_fit,";q/p (1/GeV)");
  SetPad(gPad);
  h_eQOP_fit->Draw();

  TCanvas* cPulls = new TCanvas("cPulls","",1900,800);
  cPulls->Divide(3,2,0.001,0.001);
  cPulls->cd(1);
  SetPad(gPad);
  SetHisto(h_eLOC0_fit,";d (mm)");
  h_eLOC0_fit->Draw();
  
  cPulls->cd(2);
  SetPad(gPad);
  SetHisto(h_eLOC1_fit,";z (mm)");
  h_eLOC1_fit->Draw();

  cPulls->cd(3);
  SetPad(gPad);
  SetHisto(h_eQOP_fit,";q/p (1/GeV)");
  h_eQOP_fit->Draw();
  
  cPulls->cd(4);
  SetPad(gPad);
  SetHisto(h_pull_eLOC0_fit,";pull d");
  h_pull_eLOC0_fit->Draw();

  cPulls->cd(5);
  SetPad(gPad);
  SetHisto(h_pull_eLOC1_fit,";pull z");
  h_pull_eLOC1_fit->Draw();

  cPulls->cd(6);
  SetPad(gPad);
  SetHisto(h_pull_eQOP_fit,";pull q/p");
  h_pull_eQOP_fit->Draw();
  h_pull_eQOP_fit->Fit("gaus","","",-2,2);

  cPulls->Print("pulls.png");
}