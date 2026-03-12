#include "TMath.h"

const int n = 3;
double a[n];
double b[n];
double c[n];
double sigma = 0.1;

double func(const double* xx){
  double t = xx[0];
  double alpha = xx[1];
  double d2 = 0;
  for (int i=0;i<n;i++){
    double d = a[i]*t*cos(alpha) + b[i]*t*sin(alpha) + c[i];
    d2+=d*d;
  }
  return d2;
}

double analytic(double &t, double &k, double &dt, double &dk, double& A, double& B, double& sigA2, double& sigB2, double& covAB, bool debug = 0){
  // analytical computation
  double s[n];
  double aa[n]={0};
  double bb[n]={0};
  A = 0;
  B = 0;
  sigA2 = 0;
  sigB2 = 0;
  covAB = 0;

  for (int i=0;i<n;i++){
    s[i] = sigma;
    aa[i] = 0;
    bb[i] = 0;
    for (int j=0;j<n;j++){
      double ab = (a[i]*b[j]-b[i]*a[j]);
      aa[i]+=b[j]*ab;
      bb[i]+=a[j]*ab;
    }
    A+= c[i]*aa[i];
    B+= c[i]*bb[i];
    sigA2+=s[i]*s[i]*aa[i]*aa[i];
    sigB2+=s[i]*s[i]*bb[i]*bb[i];
    covAB+=s[i]*s[i]*aa[i]*bb[i];
  }


  if (debug) printf("A=%e\n",A);
  if (debug) printf("B=%e\n",B);
  if (debug) printf("sigA2=%e\n",sigA2);  
  if (debug) printf("sigB2=%e\n",sigB2);  
  if (debug) printf("covAB=%e\n",covAB);  
  k = - B/A;
  dk = sqrt(B*B*sigA2+A*A*sigB2-2*A*B*covAB)/A/A;
  if (debug) printf("k=%e dk=%e dk/k=%f\n",k, dk, dk/k);  

  double num = 0;
  double den = 0;
  double abk[n] = {0};
  double sigT2[n] = {0};
  double sigC2[n] = {0};
  double covCT[n] = {0};
  double p1 = 0;
  double p2 = 0;
  for (int i=0;i<n;i++){
    abk[i] = a[i]+b[i]*k;
    num += c[i]*abk[i];
    den += abk[i]*abk[i];
    p1 += a[i]*b[i];
    p2 += b[i]*b[i];
  }
  t = -sqrt(1+k*k)*num/den;

  double vt[n] = {0};
  double dtdk[n] = {0};
  double dkdc[n] = {0};
  for (int i=0; i<n; i++){
    vt[i]  = -sqrt(1+k*k)*abk[i]/den;
    dtdk[i] = vt[i]*(k/(1+k*k) + b[i]/abk[i]-2*(p1+p2*k)/den);
    dkdc[i] = (aa[i]*B - bb[i]*A)/A/A;
    if (debug) printf("i=%d t[i]=%e dtdk[i]*dk=%e dtdk[i]*dk/t[i]=%f\n",i, vt[i],dtdk[i]*dk,dtdk[i]*dk/vt[i]);
  }
  double dtdc[n] = {0};
  double dt2 = 0;
  for (int i=0; i<n; i++){
    if (debug) printf("i=%d k[i]=%e dkdc[i]*dc[i]=%e dkdc[i]*dc/k[i]=%f\n",i,k,dkdc[i]*s[i],dkdc[i]*s[i]/k);
    dtdc[i] = vt[i];
    for (int j=0; j<n; j++){
      dtdc[i] += c[j]*dtdk[j]*dkdc[i];
    }
    dt2+=s[i]*s[i]*dtdc[i]*dtdc[i];
  }
  dt = sqrt(dt2);
  if (debug) printf("t=%e dt=%e dt/t=%f\n",t,dt,dt/t);

  double chi2 = 0;
  double cosa = 1./sqrt(1+k*k);
  double sina = k*cosa;
  for (int i=0;i<n;i++){
    double d = a[i]*t*cosa + b[i]*t*sina + c[i];
    chi2+=d*d/s[i]/s[i];
  }
  return chi2;
}

void analytical_space_point_maker(){
  double alphaTrueDeg = 20;
  double thetaTrue = 30*TMath::DegToRad();
  double alphaTrue = alphaTrueDeg*TMath::DegToRad();
  double tTrue = tan(thetaTrue);
  double kTrue = tan(alphaTrue);
  double r[] = { 100,  100,  100,  100,  100}; // station minimal radius
  double z[] = {2100, 2110, 2120, 2130, 2140}; // station pozition
  double phiDeg[] = {alphaTrueDeg-4, alphaTrueDeg+1, alphaTrueDeg+6, alphaTrueDeg-3, alphaTrueDeg+3};
  double d[n];
  double phi[n];
  double xTrue[n];
  double yTrue[n];
  double xmeas[n];
  double ymeas[n];
  for (int i=0;i<n;i++){
    d[i] = gRandom->Gaus(0,sigma);
    phi[i] = phiDeg[i]*TMath::DegToRad();
    xTrue[i] = z[i]*tTrue*cos(alphaTrue);
    yTrue[i] = z[i]*tTrue*sin(alphaTrue);
    xmeas[i] = xTrue[i]+d[i]*sin(phi[i]); // TODO shift to min radius
    ymeas[i] = yTrue[i]-d[i]*cos(phi[i]); // TODO shift to min radius
  }
  
  for (int i=0;i<n;i++){
    a[i] = +z[i]*sin(phi[i]);
    b[i] = -z[i]*cos(phi[i]);
    c[i] = -xmeas[i]*sin(phi[i]) + ymeas[i]*cos(phi[i]);
  }  
  
  // minimization algorithm  
  double xx[2];
  xx[0] = tTrue;
  xx[1] = alphaTrue;
  printf("%f\n",func(xx));

  ROOT::Math::Functor f(&func,2);
  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  minimum->SetMaxIterations(10000);  // for GSL
  minimum->SetTolerance(0.001);
  minimum->SetPrintLevel(1);
  minimum->SetFunction(f);
  minimum->SetVariable(0,"tan(theta)",tTrue, 0.01);
  minimum->SetVariable(1,"alpha",alphaTrue, 0.01);
  minimum->Minimize();
  auto xs = minimum->X();
  double tt = xs[0];
  double aa = xs[1];
  double kk = tan(aa);

  // analytical computation
  double t,k,dt,dk,A,B,A0,B0,sigA2,sigB2,covAB;
  analytic(t,k,dt,dk,A,B,sigA2,sigB2,covAB,1);
  
  printf("kk=%.10f k=%.10f kTrue=%f\n", kk, k, kTrue);
  printf("tt=%.10f t=%.10f tTrue=%f\n", tt, t, tTrue);
  xx[0]=t;
  xx[1]=atan(k);
  double d2 = func(xx);
  printf("Squared sum of distances: %f\n", d2);

  for (int i=0;i<n;i++){
    c[i] = -xTrue[i]*sin(phi[i]) + yTrue[i]*cos(phi[i]);
  }
  analytic(t,k,dt,dk,A0,B0,sigA2,sigB2,covAB,1);
  A0*=1e-11;
  B0*=1e-11;
  double c0true=c[0];

  // error propagation test
  TH1D* hT = new TH1D("hT","",100,tTrue-5*dt,tTrue+5*dt);
  TH1D* hK = new TH1D("hK","",100,kTrue-5*dk,kTrue+5*dk);
  TH1D* hAtanK = new TH1D("hAtanK","",100,alphaTrueDeg-0.02,alphaTrueDeg+0.02);
  TH1D* hPullK = new TH1D("hPullK","",100,-5,5);
  TH1D* hPullT = new TH1D("hPullT","",100,-5,5);
  TH2D* hKT = new TH2D("kKT",";k;t;",100,kTrue-5*dk,kTrue+5*dk,100,tTrue-5*dt,tTrue+5*dt);
  TH2D* hAB = new TH2D("hAB",";A;B;",100,A0-0.02,A0+0.02,100,B0-0.02,B0+0.02);
  TH2D* hCK = new TH2D("hCK",";C;K;",100,1-1e-2,1+1e-2, 100,1-1e-2,1+1e-2);
  TH2D* hPullKT = new TH2D("kPullKT",";pull k;pull t;",100,-5,5,100,-5,5);
  TH1D* hChi2 = new TH1D("hChi2","",1000,0,30);
  for (int ev=0;ev<1000000;ev++){
    //if (ev%100==0) printf("%d\n",ev);
    for (int i=0;i<n;i++){
      d[i] = gRandom->Gaus(0,sigma);
      xmeas[i] = xTrue[i]+d[i]*sin(phi[i]);
      ymeas[i] = yTrue[i]-d[i]*cos(phi[i]);
      c[i] = -xmeas[i]*sin(phi[i]) + ymeas[i]*cos(phi[i]);
    }  
    double chi2 = analytic(t,k,dt,dk,A,B,sigA2,sigB2,covAB);
    hT->Fill(t);
    hK->Fill(k);
    hKT->Fill(k,t);
    double pullK = (k - kTrue)/dk;
    double pullT = (t - tTrue)/dt;
    hPullK->Fill(pullK);
    hAtanK->Fill(atan(k)*TMath::RadToDeg());
    hAB->Fill(A*1.e-11,B*1.e-11);
    hCK->Fill(c[0]/c0true,k/kTrue);
    hPullT->Fill(pullT);
    hPullKT->Fill(pullK,pullT);
    hChi2->Fill(chi2);
  }
  new TCanvas;
  hK->Draw();
  hK->Fit("gaus");
  new TCanvas;
  hT->Draw();
  hT->Fit("gaus");
  new TCanvas;
  hPullK->Draw();
  hPullK->Fit("gaus");
  new TCanvas;
  hPullT->Draw();
  hPullT->Fit("gaus");
  new TCanvas;
  hAtanK->Draw();
  hAtanK->Fit("gaus");
  new TCanvas;
  hAB->Draw();
  new TCanvas;
  hCK->Draw();
  new TCanvas;
  hKT->Draw();
  new TCanvas;
  hPullKT->Draw();
  new TCanvas;
  hChi2->Draw();
}
