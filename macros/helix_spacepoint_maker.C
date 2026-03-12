#include "TMath.h"

const int n = 5;
double r[] = { 100,  100,  100,  100,  100}; // station minimal radius
double z[] = {2100, 2110, 2120, 2130, 2140}; // station pozition
double thetaTrueDeg = 30;
double alphaTrueDeg = 20;
double rTrue = 3000;
double thetaTrue = thetaTrueDeg*TMath::DegToRad();
double alphaTrue = alphaTrueDeg*TMath::DegToRad();
double phiDeg[] = {alphaTrueDeg-4., alphaTrueDeg+1., alphaTrueDeg+6., alphaTrueDeg-3., alphaTrueDeg+3.};

double d[n];
double c[n];
double s[n];
double g[n];
double sigma = 0.0001;

double helix_func(double* xx, double* par){
    double p = xx[0];
    const int nn = 5;
    double ss = 0;
    double sc = 0;
    double cc = 0;
    double gs = 0;
    double gc = 0;
    double sp[n]={0};
    double cp[n]={0};
    for (int i=0;i<nn;i++){
      double cpz = cos(p*z[i]);
      double spz = sin(p*z[i]);
      cp[i] = (c[i]- c[i]*cpz + s[i]*spz);
      sp[i] = (s[i]- c[i]*spz - s[i]*cpz);
      ss += sp[i]*sp[i];
      sc += sp[i]*cp[i];
      cc += cp[i]*cp[i];
      gs += g[i]*sp[i];
      gc += g[i]*cp[i];
    }
    double tx = (gs*sc - gc*ss);
    double ty = (gc*sc - gs*cc);
    double sum = 0;
    for (int i=0;i<nn;i++){
      sum+= (tx*cp[i] + ty*sp[i] + g[i]*(ss*cc - sc*sc))*(tx*s[i]-ty*c[i]);
    }
    return sum;
};

double helix(
  double* z, double* s, double* c, double* g, 
  double* d, double &tx, double &ty, double &p, int debug)
{
  double sp[n]={0};
  double cp[n]={0};
  auto fp = [&s,&c,&z,&g,&tx,&ty,&sp,&cp,&debug ](double p) {
    const int nn = 5;
    double ss = 0;
    double sc = 0;
    double cc = 0;
    double gs = 0;
    double gc = 0;
    for (int i=0;i<nn;i++){
      double cpz = cos(p*z[i]);
      double spz = sin(p*z[i]);
      cp[i] = (c[i]- c[i]*cpz + s[i]*spz);
      sp[i] = (s[i]- c[i]*spz - s[i]*cpz);
      if (debug>2) printf("i=%d cp=%f sp=%f\n",i, cp[i], sp[i]);
      ss += sp[i]*sp[i];
      sc += sp[i]*cp[i];
      cc += cp[i]*cp[i];
      gs += g[i]*sp[i];
      gc += g[i]*cp[i];
    }
    if (debug>1) printf("%f\n",ss);
    if (debug>1) printf("%f\n",cc);
    if (debug>1) printf("%f\n",sc);
    if (debug>1) printf("%f\n",(ss*cc - sc*sc));
    tx = (gs*sc - gc*ss);
    ty = (gc*sc - gs*cc);
    if (debug>1) printf("tx = %f\n",tx);
    if (debug>1) printf("ty = %f\n",ty);
    double sum = 0;
    for (int i=0;i<nn;i++){
      sum+= (tx*cp[i] + ty*sp[i] + g[i]*(ss*cc - sc*sc))*(tx*s[i]-ty*c[i]);
    }
    tx /=(ss*cc - sc*sc);
    ty /=(ss*cc - sc*sc);
    return sum;
  };
  if (debug) printf("fp(0)=%f\n",fp(0));
  if (debug) printf("fp(-1e-6)=%e\n",fp(-1e-6));
  if (debug) printf("fp(+1e-6)=%e\n",fp(+1e-6));

  ROOT::Math::Functor1D functor(fp);
  ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kGSL_BISECTION);
  rf.SetFunction(functor, 0.0001, 0.0002);
  rf.Solve();
  if (debug) printf("p=%e iterations=%d\n",rf.Root(),rf.Iterations());
  if (debug) printf("tx=%f\n",tx);
  if (debug) printf("ty=%f\n",ty);
  p = rf.Root();
  double chi2 = 0;
  for (int i=0;i<n;i++){
    double chi = (tx*cp[i] + ty*sp[i] + g[i])/d[i];
    chi2 += chi*chi;
  }
  if (debug) printf("chi2=%f\n",chi2);
  
  return chi2;
}

void helix_spacepoint_maker(){
  double tTrue = tan(thetaTrue);
  double pTrue = tTrue/rTrue; 
  double d[n];
  double phi[n];
  double xTrue[n];
  double yTrue[n];
  double xmeas[n];
  double ymeas[n];
  double txTrue = rTrue*cos(alphaTrue);
  double tyTrue = rTrue*sin(alphaTrue);
  printf("tTrue=%f\n",tTrue);
  printf("rTrue=%f\n",rTrue);
  printf("pTrue=%e\n",pTrue);
  printf("txTrue=%f\n",txTrue);
  printf("tyTrue=%f\n",tyTrue);  
  for (int i=0;i<n;i++){
    xTrue[i] = rTrue*(sin(pTrue*z[i])*cos(alphaTrue) + (1 - cos(pTrue*z[i]))*sin(alphaTrue));
    yTrue[i] = rTrue*(sin(pTrue*z[i])*sin(alphaTrue) - (1 - cos(pTrue*z[i]))*cos(alphaTrue));
    printf("%d %f %f %f\n",i, z[i], xTrue[i], yTrue[i]);
  }
  
  for (int i=0;i<n;i++){
    d[i] = sigma;
    double dd = gRandom->Gaus(0,d[i]);
    phi[i] = phiDeg[i]*TMath::DegToRad();
    xmeas[i] = xTrue[i] + dd*sin(phi[i]);
    ymeas[i] = yTrue[i] - dd*cos(phi[i]);
    c[i] = cos(phi[i]);
    s[i] = sin(phi[i]);
    g[i] = -xmeas[i]*sin(phi[i]) + ymeas[i]*cos(phi[i]);
  }  
  double tx, ty, p;
  double chi2 = helix(z, s, c, g, d, tx, ty, p, 1);

  TF1* f1 = new TF1("f1",&helix_func, 0.000001, 0.1, 0);
  new TCanvas;
  f1->Draw();

}
