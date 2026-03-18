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
#include "Acts/Definitions/Units.hpp"

using namespace std;

TH1D* hDeltaCotTheta2;
TH1D* hError2;
TH1D* hScatteringInRegion2;
TH1D* hMultipleScattering;
TH1D* hIm;
struct options_struct{
  float minPt = 50; // MeV
  float helixCutTolerance = 0.01;
  float sigmaScattering = 200; // 5
  float radLengthPerSeed = 0.20; // 0.20
  float bFieldInZ = 0.5 / 1000;
  float impactMax = 440;
  // derived
  float highland = 13.6 * std::sqrt(radLengthPerSeed) * (1 + 0.038 * std::log(radLengthPerSeed));
  float maxScatteringAngle = highland / minPt;
  float maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;
  float pTPerHelixRadius = Acts::UnitConstants::T * 1e6 * bFieldInZ;
  float minHelixDiameter2 = std::pow(minPt * 2 / pTPerHelixRadius, 2) * helixCutTolerance;
  float pT2perRadius = std::pow(highland / pTPerHelixRadius, 2);
  float sigmapT2perRadius = pT2perRadius * std::pow(2 * sigmaScattering, 2);
  float multipleScattering2 = maxScatteringAngle2 * std::pow(sigmaScattering, 2);
} options;


struct spacepoint{
  float sx;
  float sy;
  float sz;
  float var_r;
  float var_z;
  UInt_t event_id=1000000;
  ULong64_t geometry_id;
  ULong64_t measurement_id;

  float x() { return sx; }
  float y() { return sy; }
  float z() { return sz; }
  float radius() { return sqrt(sx*sx+sy*sy); }
  float varianceR() { return var_r; }
  float varianceZ() { return var_z; }
  int layer() { return Acts::GeometryIdentifier(geometry_id).layer(); }
};

struct LinCircle {
  LinCircle() = default;
  LinCircle(float ct, float idr, float er, float u, float v, float X, float Y)
      : cotTheta(ct), iDeltaR(idr), Er(er), U(u), V(v), x(X), y(Y) {}

  float cotTheta{0.};
  float iDeltaR{0.};
  float Er{0.};
  float U{0.};
  float V{0.};
  float x{0.};
  float y{0.};

  LinCircle(spacepoint spM, spacepoint sp, bool bottom) {
    auto xM = spM.x();
    auto yM = spM.y();
    auto zM = spM.z();    
    auto rM = spM.radius();
    auto varianceRM = spM.varianceR();
    auto varianceZM = spM.varianceZ();
    auto xSP = sp.x();
    auto ySP = sp.y();
    auto zSP = sp.z();    
    auto rSP = sp.radius();
    auto varianceRSP = sp.varianceR();
    auto varianceZSP = sp.varianceZ();
    float cosPhiM = xM / rM;
    float sinPhiM = yM / rM;
    float deltaX = xSP - xM;
    float deltaY = ySP - yM;
    float deltaZ = zSP - zM;
    float xNewFrame = deltaX * cosPhiM + deltaY * sinPhiM;
    float yNewFrame = deltaY * cosPhiM - deltaX * sinPhiM;
    float deltaR2 = (xNewFrame * xNewFrame + yNewFrame * yNewFrame);
    float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
    iDeltaR = std::sqrt(iDeltaR2);
    int bottomFactor = bottom ? -1 : 1;
    
    cotTheta = deltaZ * iDeltaR * bottomFactor;
    // conformal transformation u=x/(x²+y²) v=y/(x²+y²) transform the
    // circle into straight lines in the u/v plane the line equation can
    // be described in terms of aCoef and bCoef, where v = aCoef * u + bCoef
    U = xNewFrame * iDeltaR2;
    V = yNewFrame * iDeltaR2;
    // error term for sp-pair without correlation of middle space point
    Er = ((varianceZM + varianceZSP) + (cotTheta * cotTheta) * (varianceRM + varianceRSP)) * iDeltaR2;
    x = xNewFrame;
    y = yNewFrame;
  }

};

int filter(spacepoint spB, spacepoint spM, spacepoint spT){
  const float rM = spM.radius();
  const float cosPhiM = spM.x() / rM;
  const float sinPhiM = spM.y() / rM;
  const float varianceRM = spM.varianceR();
  const float varianceZM = spM.varianceZ();

  auto lb = LinCircle(spM,spB,1);
  float cotThetaB = lb.cotTheta;
  float Vb = lb.V;
  float Ub = lb.U;
  float ErB = lb.Er;
  float iDeltaRB = lb.iDeltaR;

  float iSinTheta2 = (1. + cotThetaB * cotThetaB);
  float sigmaSquaredPtDependent = iSinTheta2 * options.sigmapT2perRadius;
  float scatteringInRegion2 = options.multipleScattering2 * iSinTheta2;

  float sinTheta = 1 / std::sqrt(iSinTheta2);
  float cosTheta = cotThetaB * sinTheta;

  auto lt = LinCircle(spM,spT,0);
 
  float cotThetaT = lt.cotTheta;
  float rMxy = 0.;
  float ub = 0.;
  float vb = 0.;
  float ut = 0.;
  float vt = 0.;
  double rMTransf[3];
  float xB = 0.;
  float yB = 0.;
  float xT = 0.;
  float yT = 0.;
  float iDeltaRB2 = 0.;
  float iDeltaRT2 = 0.;

  float cotThetaAvg2 = cotThetaB * cotThetaT;
  if (cotThetaAvg2 <= 0) return 1;
  float error2 = lt.Er + ErB + 2 * (cotThetaAvg2 * varianceRM + varianceZM) * iDeltaRB * lt.iDeltaR;
  float deltaCotTheta = cotThetaB - cotThetaT;
  float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;
  // Fill
  // printf("%f\n",deltaCotTheta2);
  hDeltaCotTheta2->Fill(deltaCotTheta2);
  hError2->Fill(error2);
  hScatteringInRegion2->Fill(scatteringInRegion2);
  //printf("%f\n",scatteringInRegion2);

  float multipleScattering = sqrt((deltaCotTheta2)/iSinTheta2);
  printf("%f\n",multipleScattering);
  hMultipleScattering->Fill(multipleScattering);


  if (deltaCotTheta2 > (error2 + scatteringInRegion2)) return 2;

  float dU = lt.U - Ub;
  float A = (lt.V - Vb) / dU;
  float S2 = S2 = 1. + A * A;
  float B = Vb - A * Ub;
  float B2 = B * B;
  if (S2 < B2 * options.minHelixDiameter2) return 3;

  float iHelixDiameter2 = B2 / S2;
  // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta) from rad to deltaCotTheta
  float p2scatterSigma = iHelixDiameter2 * sigmaSquaredPtDependent;


  // recalculated p2scatterSigma via maxPtScattering?
  if (deltaCotTheta2 > (error2 + p2scatterSigma)) return 4;

  // A and B allow calculation of impact params in U/V plane with linear function
  // (in contrast to having to solve a quadratic function in x/y plane)
  float Im = std::abs((A - B * rM) * rM);
  hIm->Fill(Im);
  printf("Im=%f\n",Im);

  if (Im > options.impactMax) return 5;
  float curvature = B / std::sqrt(S2);
  float zOrigin = spM.z() - rM * lb.cotTheta;
  return 0;
}

void analyse_seeds(TString dir = "../build/test"){
  dir.Append("/");
  printf("%f\n",options.sigmapT2perRadius);
  
  int selectedEvent = 0;
  TFile* fSpacepoints = new TFile(dir + "spacepoints.root");
  TTree* tSpacepoints = (TTree*) fSpacepoints->Get("spacepoints");
  tSpacepoints->Print();
  spacepoint sp;
  UInt_t sevent_id;
  ULong64_t sgeometry_id;
  tSpacepoints->SetBranchAddress("x",&sp.sx);
  tSpacepoints->SetBranchAddress("y",&sp.sy);
  tSpacepoints->SetBranchAddress("z",&sp.sz);
  tSpacepoints->SetBranchAddress("var_r",&sp.var_r);
  tSpacepoints->SetBranchAddress("var_z",&sp.var_z);
  tSpacepoints->SetBranchAddress("geometry_id",&sp.geometry_id);
  tSpacepoints->SetBranchAddress("measurement_id",&sp.measurement_id);
  tSpacepoints->SetBranchAddress("event_id",&sp.event_id);

  spacepoint spB,spM,spT;
  int nAccepted = 0;
  TH1D* hRejectionReasons = new TH1D("h","h",6,0,6);
  hDeltaCotTheta2 = new TH1D("hDeltaCotTheta2","",100,0,0.2);
  hError2  = new TH1D("hError2","",100,0,0.2);
  hScatteringInRegion2  = new TH1D("hScatteringInRegion2","",1000,0,1000);
  hMultipleScattering = new TH1D("hMultipleScattering","",100,0,1);
  
  hIm = new TH1D("hIm","",100,0,400);
  
  int previous_event = -1;
  for (int is=0;is<tSpacepoints->GetEntries();is++){
    tSpacepoints->GetEntry(is);
    printf("%d %d\n",sp.layer(), sp.event_id);
    if (sp.layer()== 3) spB = sp;
    if (sp.layer()>=17 && sp.layer()<19) spM = sp;
    if (sp.layer()== 33) spT = sp;
    if (spB.event_id==previous_event) continue;
    if (spB.event_id==1000000) continue;
    if (spB.event_id!=spM.event_id) continue;
    if (spT.event_id!=spM.event_id) continue;
    printf("Event: %d\n",spB.event_id);
    previous_event = spB.event_id;
    printf("spB: x=%f y=%f z=%f\n", spB.sx, spB.sy, spB.sz);
    printf("spM: x=%f y=%f z=%f\n", spM.sx, spM.sy, spM.sz);
    printf("spT: x=%f y=%f z=%f\n", spT.sx, spT.sy, spT.sz);
    int filterResult = filter(spB,spM,spT);
    hRejectionReasons->Fill(filterResult);
    bool isAccepted = filterResult==0;
    if (isAccepted) nAccepted++;
  }
  printf("%d\n",nAccepted);
  new TCanvas;
  hDeltaCotTheta2->Draw();
  new TCanvas;
  hError2->Draw();
  new TCanvas;
  hScatteringInRegion2->Draw();
  new TCanvas;
  hMultipleScattering->Draw();
  new TCanvas;
  hIm->Draw();
  new TCanvas;
  gPad->SetLogy();
  hRejectionReasons->Draw();

}

