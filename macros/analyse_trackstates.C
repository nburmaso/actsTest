#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "vector"
#include "style.h"

  // from RootTrackStatesWriter.hpp
  enum ParameterType { ePredicted, eFiltered, eSmoothed, eUnbiased, eSize };

  /// the event number
  uint32_t m_eventNr{0};
  /// the track number
  uint32_t m_trackNr{0};

  /// Global truth hit position x
  std::vector<float>* m_t_x = new std::vector<float>;
  /// Global truth hit position y
  std::vector<float>* m_t_y = new std::vector<float>;
  /// Global truth hit position z
  std::vector<float>* m_t_z = new std::vector<float>;
  /// Global truth hit position r
  std::vector<float>* m_t_r = new std::vector<float>;
  /// Truth particle direction x at global hit position
  std::vector<float>* m_t_dx = new std::vector<float>;
  /// Truth particle direction y at global hit position
  std::vector<float>* m_t_dy = new std::vector<float>;
  /// Truth particle direction z at global hit position
  std::vector<float>* m_t_dz = new std::vector<float>;

  /// truth parameter eBoundLoc0
  std::vector<float>* m_t_eLOC0 = new std::vector<float>;
  /// truth parameter eBoundLoc1
  std::vector<float>* m_t_eLOC1 = new std::vector<float>;
  /// truth parameter ePHI
  std::vector<float>* m_t_ePHI = new std::vector<float>;
  /// truth parameter eTHETA
  std::vector<float>* m_t_eTHETA = new std::vector<float>;
  /// truth parameter eQOP
  std::vector<float>* m_t_eQOP = new std::vector<float>;
  /// truth parameter eT
  std::vector<float>* m_t_eT = new std::vector<float>;

  /// event-unique particle identifier a.k.a barcode for hits per each surface
  std::vector<std::vector<std::uint64_t>>* m_particleId = new std::vector<std::vector<std::uint64_t>>;

  /// number of all states
  unsigned int m_nStates{0};
  /// number of states with measurements
  unsigned int m_nMeasurements{0};
  /// volume identifier
  std::vector<int>* m_volumeID = new std::vector<int>;
  /// layer identifier
  std::vector<int>* m_layerID = new std::vector<int>;
  /// surface identifier
  std::vector<int>* m_moduleID = new std::vector<int>;
  /// path length
  std::vector<float>* m_pathLength = new std::vector<float>;
  /// uncalibrated measurement local x
  std::vector<float>* m_lx_hit = new std::vector<float>;
  /// uncalibrated measurement local y
  std::vector<float>* m_ly_hit = new std::vector<float>;
  /// uncalibrated measurement global x
  std::vector<float>* m_x_hit = new std::vector<float>;
  /// uncalibrated measurement global y
  std::vector<float>* m_y_hit = new std::vector<float>;
  /// uncalibrated measurement global z
  std::vector<float>* m_z_hit = new std::vector<float>;
  /// hit residual x
  std::vector<float>* m_res_x_hit = new std::vector<float>;
  /// hit residual y
  std::vector<float>* m_res_y_hit = new std::vector<float>;
  /// hit err x
  std::vector<float>* m_err_x_hit = new std::vector<float>;
  /// hit err y
  std::vector<float>* m_err_y_hit = new std::vector<float>;
  /// hit pull x
  std::vector<float>* m_pull_x_hit = new std::vector<float>;
  /// hit pull y
  std::vector<float>* m_pull_y_hit = new std::vector<float>;
  /// dimension of measurement
  std::vector<int>* m_dim_hit = new std::vector<int>;

  /// number of states which have filtered/predicted/smoothed/unbiased
  /// parameters
  std::array<int, eSize> m_nParams{};
  /// status of the filtered/predicted/smoothed/unbiased parameters
  std::array<std::vector<bool>*, eSize> m_hasParams;
  /// predicted/filtered/smoothed/unbiased parameter eLOC0
  std::array<std::vector<float>*, eSize> m_eLOC0;
  /// predicted/filtered/smoothed/unbiased parameter eLOC1
  std::array<std::vector<float>*, eSize> m_eLOC1;
  /// predicted/filtered/smoothed/unbiased parameter ePHI
  std::array<std::vector<float>*, eSize> m_ePHI;
  /// predicted/filtered/smoothed/unbiased parameter eTHETA
  std::array<std::vector<float>*, eSize> m_eTHETA;
  /// predicted/filtered/smoothed/unbiased parameter eQOP
  std::array<std::vector<float>*, eSize> m_eQOP;
  /// predicted/filtered/smoothed/unbiased parameter eT
  std::array<std::vector<float>*, eSize> m_eT;
  /// predicted/filtered/smoothed/unbiased parameter eLOC0 residual
  std::array<std::vector<float>*, eSize> m_res_eLOC0;
  /// predicted/filtered/smoothed/unbiased parameter eLOC1 residual
  std::array<std::vector<float>*, eSize> m_res_eLOC1;
  /// predicted/filtered/smoothed/unbiased parameter ePHI residual
  std::array<std::vector<float>*, eSize> m_res_ePHI;
  /// predicted/filtered/smoothed/unbiased parameter eTHETA residual
  std::array<std::vector<float>*, eSize> m_res_eTHETA;
  /// predicted/filtered/smoothed/unbiased parameter eQOP residual
  std::array<std::vector<float>*, eSize> m_res_eQOP;
  /// predicted/filtered/smoothed/unbiased parameter eT residual
  std::array<std::vector<float>*, eSize> m_res_eT;
  /// predicted/filtered/smoothed/unbiased parameter eLOC0 error
  std::array<std::vector<float>*, eSize> m_err_eLOC0;
  /// predicted/filtered/smoothed/unbiased parameter eLOC1 error
  std::array<std::vector<float>*, eSize> m_err_eLOC1;
  /// predicted/filtered/smoothed/unbiased parameter ePHI error
  std::array<std::vector<float>*, eSize> m_err_ePHI;
  /// predicted/filtered/smoothed/unbiased parameter eTHETA error
  std::array<std::vector<float>*, eSize> m_err_eTHETA;
  /// predicted/filtered/smoothed/unbiased parameter eQOP error
  std::array<std::vector<float>*, eSize> m_err_eQOP;
  /// predicted/filtered/smoothed/unbiased parameter eT error
  std::array<std::vector<float>*, eSize> m_err_eT;
  /// predicted/filtered/smoothed/unbiased parameter eLOC0 pull
  std::array<std::vector<float>*, eSize> m_pull_eLOC0;
  /// predicted/filtered/smoothed/unbiased parameter eLOC1 pull
  std::array<std::vector<float>*, eSize> m_pull_eLOC1;
  /// predicted/filtered/smoothed/unbiased parameter ePHI pull
  std::array<std::vector<float>*, eSize> m_pull_ePHI;
  /// predicted/filtered/smoothed/unbiased parameter eTHETA pull
  std::array<std::vector<float>*, eSize> m_pull_eTHETA;
  /// predicted/filtered/smoothed/unbiased parameter eQOP pull
  std::array<std::vector<float>*, eSize> m_pull_eQOP;
  /// predicted/filtered/smoothed/unbiased parameter eT pull
  std::array<std::vector<float>*, eSize> m_pull_eT;
  /// predicted/filtered/smoothed/unbiased parameter global x
  std::array<std::vector<float>*, eSize> m_x;
  /// predicted/filtered/smoothed/unbiased parameter global y
  std::array<std::vector<float>*, eSize> m_y;
  /// predicted/filtered/smoothed/unbiased parameter global z
  std::array<std::vector<float>*, eSize> m_z;
  /// predicted/filtered/smoothed/unbiased parameter px
  std::array<std::vector<float>*, eSize> m_px;
  /// predicted/filtered/smoothed/unbiased parameter py
  std::array<std::vector<float>*, eSize> m_py;
  /// predicted/filtered/smoothed/unbiased parameter pz
  std::array<std::vector<float>*, eSize> m_pz;
  /// predicted/filtered/smoothed/unbiased parameter eta
  std::array<std::vector<float>*, eSize> m_eta;
  /// predicted/filtered/smoothed/unbiased parameter pT
  std::array<std::vector<float>*, eSize> m_pT;
  /// chisq from filtering
  std::vector<float>* m_chi2 = new std::vector<float>;


void SetBranchAddresses(TTree* m_outputTree){
  for (int i=0;i<eSize;i++){
    m_hasParams[i] = new std::vector<bool>;
    m_eLOC0[i] = new std::vector<float>;
    m_eLOC1[i] = new std::vector<float>;
    m_ePHI[i] = new std::vector<float>;
    m_eTHETA[i] = new std::vector<float>;
    m_eQOP[i] = new std::vector<float>;
    m_eT[i] = new std::vector<float>;
    m_res_eLOC0[i] = new std::vector<float>;
    m_res_eLOC1[i] = new std::vector<float>;
    m_res_ePHI[i] = new std::vector<float>;
    m_res_eTHETA[i] = new std::vector<float>;
    m_res_eQOP[i] = new std::vector<float>;
    m_res_eT[i] = new std::vector<float>;
    m_err_eLOC0[i] = new std::vector<float>;
    m_err_eLOC1[i] = new std::vector<float>;
    m_err_ePHI[i] = new std::vector<float>;
    m_err_eTHETA[i] = new std::vector<float>;
    m_err_eQOP[i] = new std::vector<float>;
    m_err_eT[i] = new std::vector<float>;
    m_pull_eLOC0[i] = new std::vector<float>;
    m_pull_eLOC1[i] = new std::vector<float>;
    m_pull_ePHI[i] = new std::vector<float>;
    m_pull_eTHETA[i] = new std::vector<float>;
    m_pull_eQOP[i] = new std::vector<float>;
    m_pull_eT[i] = new std::vector<float>;
    m_x[i] = new std::vector<float>;
    m_y[i] = new std::vector<float>;
    m_z[i] = new std::vector<float>;
    m_px[i] = new std::vector<float>;
    m_py[i] = new std::vector<float>;
    m_pz[i] = new std::vector<float>;
    m_eta[i] = new std::vector<float>;
    m_pT[i] = new std::vector<float>;
  }

    // from RootTrackStatesWriter.cpp
  m_outputTree->SetBranchAddress("event_nr", &m_eventNr);
  m_outputTree->SetBranchAddress("track_nr", &m_trackNr);

  m_outputTree->SetBranchAddress("t_x", &m_t_x);
  m_outputTree->SetBranchAddress("t_y", &m_t_y);
  m_outputTree->SetBranchAddress("t_z", &m_t_z);
  m_outputTree->SetBranchAddress("t_r", &m_t_r);
  m_outputTree->SetBranchAddress("t_dx", &m_t_dx);
  m_outputTree->SetBranchAddress("t_dy", &m_t_dy);
  m_outputTree->SetBranchAddress("t_dz", &m_t_dz);
  m_outputTree->SetBranchAddress("t_eLOC0", &m_t_eLOC0);
  m_outputTree->SetBranchAddress("t_eLOC1", &m_t_eLOC1);
  m_outputTree->SetBranchAddress("t_ePHI", &m_t_ePHI);
  m_outputTree->SetBranchAddress("t_eTHETA", &m_t_eTHETA);
  m_outputTree->SetBranchAddress("t_eQOP", &m_t_eQOP);
  m_outputTree->SetBranchAddress("t_eT", &m_t_eT);
  m_outputTree->SetBranchAddress("particle_ids", &m_particleId);

  m_outputTree->SetBranchAddress("nStates", &m_nStates);
  m_outputTree->SetBranchAddress("nMeasurements", &m_nMeasurements);
  m_outputTree->SetBranchAddress("volume_id", &m_volumeID);
  m_outputTree->SetBranchAddress("layer_id", &m_layerID);
  m_outputTree->SetBranchAddress("module_id", &m_moduleID);
  m_outputTree->SetBranchAddress("pathLength", &m_pathLength);
  m_outputTree->SetBranchAddress("l_x_hit", &m_lx_hit);
  m_outputTree->SetBranchAddress("l_y_hit", &m_ly_hit);
  m_outputTree->SetBranchAddress("g_x_hit", &m_x_hit);
  m_outputTree->SetBranchAddress("g_y_hit", &m_y_hit);
  m_outputTree->SetBranchAddress("g_z_hit", &m_z_hit);
  m_outputTree->SetBranchAddress("res_x_hit", &m_res_x_hit);
  m_outputTree->SetBranchAddress("res_y_hit", &m_res_y_hit);
  m_outputTree->SetBranchAddress("err_x_hit", &m_err_x_hit);
  m_outputTree->SetBranchAddress("err_y_hit", &m_err_y_hit);
  m_outputTree->SetBranchAddress("pull_x_hit", &m_pull_x_hit);
  m_outputTree->SetBranchAddress("pull_y_hit", &m_pull_y_hit);
  m_outputTree->SetBranchAddress("dim_hit", &m_dim_hit);

  m_outputTree->SetBranchAddress("nPredicted", &m_nParams[ePredicted]);
  m_outputTree->SetBranchAddress("predicted", &m_hasParams[ePredicted]);
  m_outputTree->SetBranchAddress("eLOC0_prt", &m_eLOC0[ePredicted]);
  m_outputTree->SetBranchAddress("eLOC1_prt", &m_eLOC1[ePredicted]);
  m_outputTree->SetBranchAddress("ePHI_prt", &m_ePHI[ePredicted]);
  m_outputTree->SetBranchAddress("eTHETA_prt", &m_eTHETA[ePredicted]);
  m_outputTree->SetBranchAddress("eQOP_prt", &m_eQOP[ePredicted]);
  m_outputTree->SetBranchAddress("eT_prt", &m_eT[ePredicted]);
  m_outputTree->SetBranchAddress("res_eLOC0_prt", &m_res_eLOC0[ePredicted]);
  m_outputTree->SetBranchAddress("res_eLOC1_prt", &m_res_eLOC1[ePredicted]);
  m_outputTree->SetBranchAddress("res_ePHI_prt", &m_res_ePHI[ePredicted]);
  m_outputTree->SetBranchAddress("res_eTHETA_prt", &m_res_eTHETA[ePredicted]);
  m_outputTree->SetBranchAddress("res_eQOP_prt", &m_res_eQOP[ePredicted]);
  m_outputTree->SetBranchAddress("res_eT_prt", &m_res_eT[ePredicted]);
  m_outputTree->SetBranchAddress("err_eLOC0_prt", &m_err_eLOC0[ePredicted]);
  m_outputTree->SetBranchAddress("err_eLOC1_prt", &m_err_eLOC1[ePredicted]);
  m_outputTree->SetBranchAddress("err_ePHI_prt", &m_err_ePHI[ePredicted]);
  m_outputTree->SetBranchAddress("err_eTHETA_prt", &m_err_eTHETA[ePredicted]);
  m_outputTree->SetBranchAddress("err_eQOP_prt", &m_err_eQOP[ePredicted]);
  m_outputTree->SetBranchAddress("err_eT_prt", &m_err_eT[ePredicted]);
  m_outputTree->SetBranchAddress("pull_eLOC0_prt", &m_pull_eLOC0[ePredicted]);
  m_outputTree->SetBranchAddress("pull_eLOC1_prt", &m_pull_eLOC1[ePredicted]);
  m_outputTree->SetBranchAddress("pull_ePHI_prt", &m_pull_ePHI[ePredicted]);
  m_outputTree->SetBranchAddress("pull_eTHETA_prt", &m_pull_eTHETA[ePredicted]);
  m_outputTree->SetBranchAddress("pull_eQOP_prt", &m_pull_eQOP[ePredicted]);
  m_outputTree->SetBranchAddress("pull_eT_prt", &m_pull_eT[ePredicted]);
  m_outputTree->SetBranchAddress("g_x_prt", &m_x[ePredicted]);
  m_outputTree->SetBranchAddress("g_y_prt", &m_y[ePredicted]);
  m_outputTree->SetBranchAddress("g_z_prt", &m_z[ePredicted]);
  m_outputTree->SetBranchAddress("px_prt", &m_px[ePredicted]);
  m_outputTree->SetBranchAddress("py_prt", &m_py[ePredicted]);
  m_outputTree->SetBranchAddress("pz_prt", &m_pz[ePredicted]);
  m_outputTree->SetBranchAddress("eta_prt", &m_eta[ePredicted]);
  m_outputTree->SetBranchAddress("pT_prt", &m_pT[ePredicted]);

  m_outputTree->SetBranchAddress("nFiltered", &m_nParams[eFiltered]);
  m_outputTree->SetBranchAddress("filtered", &m_hasParams[eFiltered]);
  m_outputTree->SetBranchAddress("eLOC0_flt", &m_eLOC0[eFiltered]);
  m_outputTree->SetBranchAddress("eLOC1_flt", &m_eLOC1[eFiltered]);
  m_outputTree->SetBranchAddress("ePHI_flt", &m_ePHI[eFiltered]);
  m_outputTree->SetBranchAddress("eTHETA_flt", &m_eTHETA[eFiltered]);
  m_outputTree->SetBranchAddress("eQOP_flt", &m_eQOP[eFiltered]);
  m_outputTree->SetBranchAddress("eT_flt", &m_eT[eFiltered]);
  m_outputTree->SetBranchAddress("res_eLOC0_flt", &m_res_eLOC0[eFiltered]);
  m_outputTree->SetBranchAddress("res_eLOC1_flt", &m_res_eLOC1[eFiltered]);
  m_outputTree->SetBranchAddress("res_ePHI_flt", &m_res_ePHI[eFiltered]);
  m_outputTree->SetBranchAddress("res_eTHETA_flt", &m_res_eTHETA[eFiltered]);
  m_outputTree->SetBranchAddress("res_eQOP_flt", &m_res_eQOP[eFiltered]);
  m_outputTree->SetBranchAddress("res_eT_flt", &m_res_eT[eFiltered]);
  m_outputTree->SetBranchAddress("err_eLOC0_flt", &m_err_eLOC0[eFiltered]);
  m_outputTree->SetBranchAddress("err_eLOC1_flt", &m_err_eLOC1[eFiltered]);
  m_outputTree->SetBranchAddress("err_ePHI_flt", &m_err_ePHI[eFiltered]);
  m_outputTree->SetBranchAddress("err_eTHETA_flt", &m_err_eTHETA[eFiltered]);
  m_outputTree->SetBranchAddress("err_eQOP_flt", &m_err_eQOP[eFiltered]);
  m_outputTree->SetBranchAddress("err_eT_flt", &m_err_eT[eFiltered]);
  m_outputTree->SetBranchAddress("pull_eLOC0_flt", &m_pull_eLOC0[eFiltered]);
  m_outputTree->SetBranchAddress("pull_eLOC1_flt", &m_pull_eLOC1[eFiltered]);
  m_outputTree->SetBranchAddress("pull_ePHI_flt", &m_pull_ePHI[eFiltered]);
  m_outputTree->SetBranchAddress("pull_eTHETA_flt", &m_pull_eTHETA[eFiltered]);
  m_outputTree->SetBranchAddress("pull_eQOP_flt", &m_pull_eQOP[eFiltered]);
  m_outputTree->SetBranchAddress("pull_eT_flt", &m_pull_eT[eFiltered]);
  m_outputTree->SetBranchAddress("g_x_flt", &m_x[eFiltered]);
  m_outputTree->SetBranchAddress("g_y_flt", &m_y[eFiltered]);
  m_outputTree->SetBranchAddress("g_z_flt", &m_z[eFiltered]);
  m_outputTree->SetBranchAddress("px_flt", &m_px[eFiltered]);
  m_outputTree->SetBranchAddress("py_flt", &m_py[eFiltered]);
  m_outputTree->SetBranchAddress("pz_flt", &m_pz[eFiltered]);
  m_outputTree->SetBranchAddress("eta_flt", &m_eta[eFiltered]);
  m_outputTree->SetBranchAddress("pT_flt", &m_pT[eFiltered]);

  m_outputTree->SetBranchAddress("nSmoothed", &m_nParams[eSmoothed]);
  m_outputTree->SetBranchAddress("smoothed", &m_hasParams[eSmoothed]);
  m_outputTree->SetBranchAddress("eLOC0_smt", &m_eLOC0[eSmoothed]);
  m_outputTree->SetBranchAddress("eLOC1_smt", &m_eLOC1[eSmoothed]);
  m_outputTree->SetBranchAddress("ePHI_smt", &m_ePHI[eSmoothed]);
  m_outputTree->SetBranchAddress("eTHETA_smt", &m_eTHETA[eSmoothed]);
  m_outputTree->SetBranchAddress("eQOP_smt", &m_eQOP[eSmoothed]);
  m_outputTree->SetBranchAddress("eT_smt", &m_eT[eSmoothed]);
  m_outputTree->SetBranchAddress("res_eLOC0_smt", &m_res_eLOC0[eSmoothed]);
  m_outputTree->SetBranchAddress("res_eLOC1_smt", &m_res_eLOC1[eSmoothed]);
  m_outputTree->SetBranchAddress("res_ePHI_smt", &m_res_ePHI[eSmoothed]);
  m_outputTree->SetBranchAddress("res_eTHETA_smt", &m_res_eTHETA[eSmoothed]);
  m_outputTree->SetBranchAddress("res_eQOP_smt", &m_res_eQOP[eSmoothed]);
  m_outputTree->SetBranchAddress("res_eT_smt", &m_res_eT[eSmoothed]);
  m_outputTree->SetBranchAddress("err_eLOC0_smt", &m_err_eLOC0[eSmoothed]);
  m_outputTree->SetBranchAddress("err_eLOC1_smt", &m_err_eLOC1[eSmoothed]);
  m_outputTree->SetBranchAddress("err_ePHI_smt", &m_err_ePHI[eSmoothed]);
  m_outputTree->SetBranchAddress("err_eTHETA_smt", &m_err_eTHETA[eSmoothed]);
  m_outputTree->SetBranchAddress("err_eQOP_smt", &m_err_eQOP[eSmoothed]);
  m_outputTree->SetBranchAddress("err_eT_smt", &m_err_eT[eSmoothed]);
  m_outputTree->SetBranchAddress("pull_eLOC0_smt", &m_pull_eLOC0[eSmoothed]);
  m_outputTree->SetBranchAddress("pull_eLOC1_smt", &m_pull_eLOC1[eSmoothed]);
  m_outputTree->SetBranchAddress("pull_ePHI_smt", &m_pull_ePHI[eSmoothed]);
  m_outputTree->SetBranchAddress("pull_eTHETA_smt", &m_pull_eTHETA[eSmoothed]);
  m_outputTree->SetBranchAddress("pull_eQOP_smt", &m_pull_eQOP[eSmoothed]);
  m_outputTree->SetBranchAddress("pull_eT_smt", &m_pull_eT[eSmoothed]);
  m_outputTree->SetBranchAddress("g_x_smt", &m_x[eSmoothed]);
  m_outputTree->SetBranchAddress("g_y_smt", &m_y[eSmoothed]);
  m_outputTree->SetBranchAddress("g_z_smt", &m_z[eSmoothed]);
  m_outputTree->SetBranchAddress("px_smt", &m_px[eSmoothed]);
  m_outputTree->SetBranchAddress("py_smt", &m_py[eSmoothed]);
  m_outputTree->SetBranchAddress("pz_smt", &m_pz[eSmoothed]);
  m_outputTree->SetBranchAddress("eta_smt", &m_eta[eSmoothed]);
  m_outputTree->SetBranchAddress("pT_smt", &m_pT[eSmoothed]);

  m_outputTree->SetBranchAddress("nUnbiased", &m_nParams[eUnbiased]);
  m_outputTree->SetBranchAddress("unbiased", &m_hasParams[eUnbiased]);
  m_outputTree->SetBranchAddress("eLOC0_ubs", &m_eLOC0[eUnbiased]);
  m_outputTree->SetBranchAddress("eLOC1_ubs", &m_eLOC1[eUnbiased]);
  m_outputTree->SetBranchAddress("ePHI_ubs", &m_ePHI[eUnbiased]);
  m_outputTree->SetBranchAddress("eTHETA_ubs", &m_eTHETA[eUnbiased]);
  m_outputTree->SetBranchAddress("eQOP_ubs", &m_eQOP[eUnbiased]);
  m_outputTree->SetBranchAddress("eT_ubs", &m_eT[eUnbiased]);
  m_outputTree->SetBranchAddress("res_eLOC0_ubs", &m_res_eLOC0[eUnbiased]);
  m_outputTree->SetBranchAddress("res_eLOC1_ubs", &m_res_eLOC1[eUnbiased]);
  m_outputTree->SetBranchAddress("res_ePHI_ubs", &m_res_ePHI[eUnbiased]);
  m_outputTree->SetBranchAddress("res_eTHETA_ubs", &m_res_eTHETA[eUnbiased]);
  m_outputTree->SetBranchAddress("res_eQOP_ubs", &m_res_eQOP[eUnbiased]);
  m_outputTree->SetBranchAddress("res_eT_ubs", &m_res_eT[eUnbiased]);
  m_outputTree->SetBranchAddress("err_eLOC0_ubs", &m_err_eLOC0[eUnbiased]);
  m_outputTree->SetBranchAddress("err_eLOC1_ubs", &m_err_eLOC1[eUnbiased]);
  m_outputTree->SetBranchAddress("err_ePHI_ubs", &m_err_ePHI[eUnbiased]);
  m_outputTree->SetBranchAddress("err_eTHETA_ubs", &m_err_eTHETA[eUnbiased]);
  m_outputTree->SetBranchAddress("err_eQOP_ubs", &m_err_eQOP[eUnbiased]);
  m_outputTree->SetBranchAddress("err_eT_ubs", &m_err_eT[eUnbiased]);
  m_outputTree->SetBranchAddress("pull_eLOC0_ubs", &m_pull_eLOC0[eUnbiased]);
  m_outputTree->SetBranchAddress("pull_eLOC1_ubs", &m_pull_eLOC1[eUnbiased]);
  m_outputTree->SetBranchAddress("pull_ePHI_ubs", &m_pull_ePHI[eUnbiased]);
  m_outputTree->SetBranchAddress("pull_eTHETA_ubs", &m_pull_eTHETA[eUnbiased]);
  m_outputTree->SetBranchAddress("pull_eQOP_ubs", &m_pull_eQOP[eUnbiased]);
  m_outputTree->SetBranchAddress("pull_eT_ubs", &m_pull_eT[eUnbiased]);
  m_outputTree->SetBranchAddress("g_x_ubs", &m_x[eUnbiased]);
  m_outputTree->SetBranchAddress("g_y_ubs", &m_y[eUnbiased]);
  m_outputTree->SetBranchAddress("g_z_ubs", &m_z[eUnbiased]);
  m_outputTree->SetBranchAddress("px_ubs", &m_px[eUnbiased]);
  m_outputTree->SetBranchAddress("py_ubs", &m_py[eUnbiased]);
  m_outputTree->SetBranchAddress("pz_ubs", &m_pz[eUnbiased]);
  m_outputTree->SetBranchAddress("eta_ubs", &m_eta[eUnbiased]);
  m_outputTree->SetBranchAddress("pT_ubs", &m_pT[eUnbiased]);

  m_outputTree->SetBranchAddress("chi2", &m_chi2);
}

void analyse_trackstates(){
  TFile* f = new TFile("build/trackstates_ckf.root");
  TTree* t = (TTree*) f->Get("trackstates");
  SetBranchAddresses(t);

  TH1D* h_eQOP_fit = new TH1D("h_eQOP_fit","",200,0,2);
  TH1D* h_eLOC0_fit = new TH1D("h_eLOC0_fit","",200,-90,90);
  TH1D* h_eLOC1_fit = new TH1D("h_eLOC1_fit","",200,-90,90);
  TH1D* h_pull_eQOP_fit = new TH1D("h_pull_eQOP_fit","",200,-10,10);
  TH1D* h_pull_eLOC0_fit = new TH1D("h_pull_eLOC0_fit","",200,-10,10);
  TH1D* h_pull_eLOC1_fit = new TH1D("h_pull_eLOC1_fit","",200,-10,10);

  TH1D* h_pT_ubs = new TH1D("h_pT_ubs","",200,0,0.4);

  int nEvents = t->GetEntries();
  for (int ev=0; ev<nEvents; ev++){
    t->GetEntry(ev);
    // fill unbiased parameters and pulls at the first measurement
    int im = 0;
    if (m_hasParams[eUnbiased]->at(im)) {
      h_eQOP_fit->Fill(m_eQOP[eUnbiased]->at(im));
      h_eLOC0_fit->Fill(m_res_eLOC0[eUnbiased]->at(im));
      h_eLOC1_fit->Fill(m_res_eLOC1[eUnbiased]->at(im));
      h_pull_eQOP_fit->Fill(m_pull_eQOP[eUnbiased]->at(im));
      h_pull_eLOC0_fit->Fill(m_pull_eLOC0[eUnbiased]->at(im));
      h_pull_eLOC1_fit->Fill(m_pull_eLOC1[eUnbiased]->at(im));
      h_pT_ubs->Fill(m_pT[eUnbiased]->at(im));
    }
  }

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

  cPulls->Print("pulls_mes0_unbiased.png");

  new TCanvas;
  h_pT_ubs->Draw();
}