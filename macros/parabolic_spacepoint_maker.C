#include "TMath.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <complex>

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
double kk = 0;
double tx = 0;
double ty = 0;




const double EPS = 1e-12;

std::vector<double> solveCubic(double a, double b, double c, double d) {
    std::vector<double> roots;
    if (std::abs(a) < EPS) { // не кубическое, а квадратное
        if (std::abs(b) < EPS) { 
            if (std::abs(c) < EPS) return roots; // нет уравнения
            roots.push_back(-d/c);
            return roots;
        }
        double D = c*c - 4*b*d;
        if (D >= 0) {
            roots.push_back((-c + std::sqrt(D))/(2*b));
            roots.push_back((-c - std::sqrt(D))/(2*b));
        }
        return roots;
    }

    // нормализация
    double p = b/a;
    double q = c/a;
    double r = d/a;

    // приводим к форме x^3 + px^2 + qx + r = 0
    double a1 = q - b*b/(3*a*a);
    double a0 = 2*b*b*b/(27*a*a*a) - b*c/(3*a*a) + d/a;

    double Q = a1/3.0;
    double R = -a0/2.0;
    double D = R*R + Q*Q*Q; // дискриминант

    if (D > EPS) {
        // один реальный корень
        double S = std::cbrt(R + std::sqrt(D));
        double T = std::cbrt(R - std::sqrt(D));
        double x = S + T - b/(3*a);
        roots.push_back(x);
    } else if (std::abs(D) <= EPS) {
        // три корня, хотя бы два совпадают
        double S = std::cbrt(R);
        roots.push_back(2*S - b/(3*a));
        roots.push_back(-S - b/(3*a));
    } else {
        // три различных вещественных корня
        double theta = std::acos(R / std::sqrt(-Q*Q*Q));
        double sqrtQ = 2*std::sqrt(-Q);
        roots.push_back(sqrtQ * std::cos(theta/3.0) - b/(3*a));
        roots.push_back(sqrtQ * std::cos((theta + 2*M_PI)/3.0) - b/(3*a));
        roots.push_back(sqrtQ * std::cos((theta + 4*M_PI)/3.0) - b/(3*a));
    }

    return roots;
}

// // Solve cubic using Cardano formula
// std::vector<double> solveCubic(double a, double b, double c, double d) {
//     std::vector<double> roots;

//     if (std::abs(a) < 1e-12) { // degenerate to quadratic
//         double D = c*c - 4*b*d;
//         if (D >= 0) {
//             roots.push_back((-c + std::sqrt(D)) / (2*b));
//             roots.push_back((-c - std::sqrt(D)) / (2*b));
//         }
//         return roots;
//     }

//     // normalize cubic: x^3 + p x^2 + q x + r = 0
//     double p = b / a;
//     double q = c / a;
//     double r = d / a;

//     double p3 = p / 3.0;
//     double Q = (3*q - p*p) / 9.0;
//     double R = (9*p*q - 27*r - 2*p*p*p) / 54.0;
//     double D = Q*Q*Q + R*R; // discriminant

//     if (D >= 0) { // one real root
//         double S = std::cbrt(R + std::sqrt(D));
//         double T = std::cbrt(R - std::sqrt(D));
//         double root = -p3 + S + T;
//         roots.push_back(root);
//     } else { // three real roots
//         double theta = std::acos(R / std::sqrt(-Q*Q*Q));
//         double sqrtQ = std::sqrt(-Q);
//         roots.push_back(2*sqrtQ*std::cos(theta/3) - p3);
//         roots.push_back(2*sqrtQ*std::cos((theta+2*M_PI)/3) - p3);
//         roots.push_back(2*sqrtQ*std::cos((theta+4*M_PI)/3) - p3);
//     }

//     return roots;
// }

struct Measurement {
    double z;
    double phi;
    double g;
};

struct TXTYK {
    double tx;
    double ty;
    double k;
};

// Compute sums needed for cubic
TXTYK fitHelixSecondOrder(const std::vector<Measurement>& meas) {
    int n = meas.size();

    std::vector<double> A(n), B(n), C(n), D(n);
    std::vector<double> s(n), c(n), z(n), g(n);

    for(int i=0;i<n;i++){
        z[i] = meas[i].z;
        s[i] = std::sin(meas[i].phi);
        c[i] = std::cos(meas[i].phi);
        g[i] = meas[i].g;
        printf("%d %f %f %f %f\n",i, z[i],s[i],c[i],g[i]);

        A[i] = z[i]*s[i];
        B[i] = -z[i]*c[i];
        C[i] = -z[i]*B[i];
        D[i] = z[i]*A[i];
    }
    // tx=3.766159e-01;
    // ty=-2.444310e-01;
    // kk=-4.203342e-04;

    double sum_tx = 0;
    double sum_ty = 0;
    printf("check %f %f %f\n",tx, ty, kk);
    for (int i=0;i<n;i++){
      double a = A[i] + kk*C[i];
      double b = B[i] + kk*D[i];
      double d = a*tx + b*ty + g[i];
      sum_tx += d*a;
      sum_ty += d*b;
    }
    printf("sum_tx=%f\n",sum_tx);
    printf("sum_ty=%f\n",sum_ty);

    // Precompute sums for cubic coefficients
    double sumAA0=0, sumAA1=0, sumAA2=0;
    double sumBB0=0, sumBB1=0, sumBB2=0;
    double sumAB0=0, sumAB1=0, sumAB2=0;
    double sumGA0=0, sumGA1=0;
    double sumGB0=0, sumGB1=0;

    for(int i=0;i<n;i++){
        sumAA0 += A[i]*A[i]; sumAA1 += A[i]*C[i];           sumAA2 += C[i]*C[i];
        sumBB0 += B[i]*B[i]; sumBB1 += B[i]*D[i];           sumBB2 += D[i]*D[i];
        sumAB0 += A[i]*B[i]; sumAB1 += A[i]*D[i]+B[i]*C[i]; sumAB2 += C[i]*D[i];
        sumGA0 += g[i]*A[i]; sumGA1 += g[i]*C[i];
        sumGB0 += g[i]*B[i]; sumGB1 += g[i]*D[i];
    }

    // Cubic coefficients for k: alpha3 k^3 + alpha2 k^2 + alpha1 k + alpha0 = 0
    double alpha3 = sumAA2*sumBB1 - sumBB2*sumAB1 - sumAB2*sumAB1;
    double alpha2 = sumAA2*sumBB0 - sumBB2*sumAA0 + sumAA1*sumBB1 - sumAB1*sumAB1 - sumGA1*sumAB1 + sumGB1*sumAB1;
    double alpha1 = sumAA1*sumBB0 - sumBB1*sumAA0 - sumAB0*sumAB1 - sumGA0*sumAB1 - sumGA1*sumAB0 + sumGB0*sumAB1 + sumGB1*sumAB0;
    double alpha0 = sumAA0*sumBB0 - sumAB0*sumAB0 - sumGA0*sumAB0 + sumGB0*sumAB0;

    // Solve cubic
    std::vector<double> k_roots = solveCubic(alpha3, alpha2, alpha1, alpha0);
    double k0 = k_roots[0];
    printf("%f %f\n", k0, alpha3*k0*k0*k0+alpha2*k0*k0+ alpha1*k0 + alpha0);
    k0 = k_roots[1];
    printf("%f %f\n", k0, alpha3*k0*k0*k0+alpha2*k0*k0+ alpha1*k0 + alpha0);
    k0 = k_roots[2];
    printf("%f %f\n", k0, alpha3*k0*k0*k0+alpha2*k0*k0+ alpha1*k0 + alpha0);
    printf("alpha0 = %f\n", alpha0);

    // Pick the physically meaningful root (e.g., smallest abs)
    double k = *std::min_element(k_roots.begin(), k_roots.end(), [](double a,double b){return std::abs(a)<std::abs(b);});

    // Compute t_x, t_y for this k
    double AA = sumAA0 + 2*k*sumAA1 + k*k*sumAA2;
    double BB = sumBB0 + 2*k*sumBB1 + k*k*sumBB2;
    double AB = sumAB0 + k*sumAB1 + k*k*sumAB2;
    double GA = sumGA0 + k*sumGA1;
    double GB = sumGB0 + k*sumGB1;

    double delta = AA*BB - AB*AB;

    double res_tx = (GB*AB - GA*BB)/delta;
    double res_ty = (GA*AB - GB*AA)/delta;

    return {res_tx, res_ty, k};
}


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
double sigma = 0.001;

double func(const double* xx){
  double k = xx[0];
  double tx = xx[1];
  double ty = xx[2];
  double d2 = 0;
  for (int i=0;i<n;i++){
    double sk = s[i] + k*z[i]*c[i];
    double ck = c[i] - k*z[i]*s[i];
    double d = sk*tx - ck*ty + g[i];
    // printf("%f %f %f %f %f %f\n",d,tx,ty,sk,ck,g[i]);
    d2+=d*d;
  }
  return d2;
}

double func2(double* xx, double* par){
  xx[1] = par[0];
  xx[2] = par[1];
  return func(xx);
}

double func1(double* xx, double* par){
  double k = xx[0];
  double sk[n]={0};
  double ck[n]={0};  
  double ss = 0;
  double sc = 0;
  double cc = 0;
  double gs = 0;
  double gc = 0;
  for (int i=0;i<n;i++){
    sk[i] = s[i] + k*z[i]*c[i];
    ck[i] = c[i] - k*z[i]*s[i];
    ss += sk[i]*sk[i];
    sc += sk[i]*ck[i];
    cc += ck[i]*ck[i];
    gs += g[i]*sk[i];
    gc += g[i]*ck[i];
  }
  double tx = (gc*sc - gs*cc);
  double ty = (gc*ss - gs*sc);

  double d = 0;
  for (int i=0;i<n;i++){
    d+=(tx*sk[i]-ty*ck[i]+g[i]*(ss*cc - sc*sc))*(ty*z[i]*s[i]+tx*z[i]*c[i]);
  }

  tx/=(ss*cc - sc*sc);
  ty/=(ss*cc - sc*sc);
  return d;
}

double func_k(double k){
  double sk[n]={0};
  double ck[n]={0};  
  double ss = 0;
  double sc = 0;
  double cc = 0;
  double gs = 0;
  double gc = 0;
  for (int i=0;i<n;i++){
    sk[i] = s[i] + k*z[i]*c[i];
    ck[i] = c[i] - k*z[i]*s[i];
    ss += sk[i]*sk[i];
    sc += sk[i]*ck[i];
    cc += ck[i]*ck[i];
    gs += g[i]*sk[i];
    gc += g[i]*ck[i];
  }
  double tx = (gc*sc - gs*cc)/(ss*cc - sc*sc);
  double ty = (gc*ss - gs*sc)/(ss*cc - sc*sc);

  double d = 0;
  for (int i=0;i<n;i++){
    d+=(tx*sk[i]-ty*ck[i]+g[i])*(ty*z[i]*s[i]+tx*z[i]*c[i]);
  }

  return d;
}


double parabola(
  double* z, double* s, double* c, double* g, 
  double* d, double &tx, double &ty, double &p, bool debug)
{
  double sk[n]={0};
  double ck[n]={0};  
  auto fk = [&s,&c,&z,&g,&tx,&ty,&sk,&ck,&debug ](double k) {
    const int nn = 5;
    double ss = 0;
    double sc = 0;
    double cc = 0;
    double gs = 0;
    double gc = 0;
    for (int i=0;i<nn;i++){
      sk[i] = (s[i] + k*z[i]*c[i]);
      ck[i] = (c[i] - k*z[i]*s[i]);
      ss += sk[i]*sk[i];
      sc += sk[i]*ck[i];
      cc += ck[i]*ck[i];
      gs += g[i]*sk[i];
      gc += g[i]*ck[i];
    }
    tx = (gc*sc - gs*cc)/(ss*cc - sc*sc);
    ty = (gc*ss - gs*sc)/(ss*cc - sc*sc);
    if (debug) printf("tx = %f\n",tx);
    if (debug) printf("ty = %f\n",ty);
    double sum = 0;
    for (int i=0;i<nn;i++){
      sum+=(tx*sk[i]-ty*ck[i]+g[i])*(ty*z[i]*s[i]+tx*z[i]*c[i]);
    }
    return sum;
  };
  // printf("fk(0)=%f\n",fk(0));

  ROOT::Math::Functor1D functor(fk);
  ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kGSL_BISECTION);
  rf.SetFunction(functor, -0.0002, 0.0002);
  rf.Solve();
  printf("k=%e iterations=%d\n",rf.Root(),rf.Iterations());

  p = 2*rf.Root();
  double chi2 = 0;
  for (int i=0;i<n;i++){
    double chi = (sk[i]*tx - ck[i]*ty + g[i])/d[i];
    chi2 += chi*chi;
  }

  double sum_tx = 0;
  double sum_ty = 0;
  kk = rf.Root();
  printf("fk(kk)=%f\n",fk(kk));
  printf("check %f %f %f\n",tx, ty, kk);
  for (int i=0;i<n;i++){
    double a = s[i] + kk*z[i]*c[i];
    double b = c[i] - kk*z[i]*s[i];
    printf("a=%e sk=%e\n",a,sk[i]);
    // a = sk[i];
    // b = ck[i];
    double d = a*tx - b*ty + g[i];
    sum_tx += d*a;
    sum_ty += d*b;
  }
  printf("sum_tx=%f\n",sum_tx);
  printf("sum_ty=%f\n",sum_ty);


  return chi2;
}


void parabolic_spacepoint_maker(){
  double tTrue = tan(thetaTrue);
  double pTrue = tTrue/rTrue; 
  double d[n];
  double phi[n];
  double xTrue[n];
  double yTrue[n];
  double xmeas[n];
  double ymeas[n];
  double txTrue = rTrue*pTrue*cos(alphaTrue);
  double tyTrue = rTrue*pTrue*sin(alphaTrue);
  double kTrue = pTrue/2.;
  printf("tTrue=%f\n",tTrue);
  printf("rTrue=%f\n",rTrue);
  printf("pTrue=%f\n",pTrue);
  printf("txTrue=%f\n",txTrue);
  printf("tyTrue=%f\n",tyTrue);  
  printf("kTrue=%e\n",kTrue);  
  for (int i=0;i<n;i++){
    xTrue[i] = txTrue*z[i] + kTrue*tyTrue*z[i]*z[i];
    yTrue[i] = tyTrue*z[i] - kTrue*txTrue*z[i]*z[i];
    printf("%d %f %f %f\n",i, z[i], xTrue[i], yTrue[i]);
  }
  
  for (int i=0;i<n;i++){
    d[i] = sigma;
    double dd = gRandom->Gaus(0,d[i]);
    phi[i] = phiDeg[i]*TMath::DegToRad();
    xmeas[i] = xTrue[i] + dd*sin(phi[i]);
    ymeas[i] = yTrue[i] - dd*cos(phi[i]);
    c[i] = z[i]*cos(phi[i]);
    s[i] = z[i]*sin(phi[i]);
    g[i] = -xmeas[i]*sin(phi[i]) + ymeas[i]*cos(phi[i]);
    printf("%d %f %f %f %f\n",i, z[i],sin(phi[i]),cos(phi[i]),g[i]);
  }  

  
  // minimization algorithm  
  double xx[3];
  xx[0] = kTrue;
  xx[1] = txTrue;
  xx[2] = tyTrue;
  printf("%f\n",func(xx));

  ROOT::Math::Functor f(&func,3);
  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  minimum->SetMaxIterations(10000);  // for GSL
  minimum->SetTolerance(0.0001);
  minimum->SetPrintLevel(1);
  minimum->SetFunction(f);
  minimum->SetVariable(0,"k",kTrue, 0.01);
  minimum->SetVariable(1,"tx",txTrue, 0.01);
  minimum->SetVariable(2,"ty",tyTrue, 0.01);
  // minimum->FixVariable(1);
  // minimum->FixVariable(2);
  minimum->Minimize();
  auto xs = minimum->X();
  kk = xs[0];
  tx = xs[1];
  ty = xs[2];
  double alpha = atan(ty/tx);
  double t = sqrt(tx*tx+ty*ty);
  double r = t/2./kk;
  printf("alpha=%f\n",alpha*TMath::RadToDeg());
  printf("t=%f\n",t);
  printf("theta=%f\n",atan(t)*TMath::RadToDeg());
  printf("r=%f\n",r);


  new TCanvas;
  TF1* f2 = new TF1("f2",&func2,0.,0.0002,2);
  f2->SetParameter(0,tx);
  f2->SetParameter(1,ty);
  f2->Draw();
  TH1* h = f2->CreateHistogram();
  h->Draw();
  h->Fit("pol2");

  new TCanvas;
  TF1* f1 = new TF1("f1",&func1,-0.02,0.02,0);
  f1->SetNpx(1000);
  f1->Draw();


  ROOT::Math::Functor1D fk(&func_k);
  ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kGSL_BISECTION);
  rf.SetFunction(fk, -0.0002, 0.0002);
  rf.Solve();

  printf("k=%e iterations=%d\n",rf.Root(),rf.Iterations());

  double p;
  double chi2 = parabola(z, s, c, g, d, tx, ty, p, 0);

  std::vector<Measurement> meas;
  for (int i=0;i<n;i++){
    meas.emplace_back(z[i],phi[i],g[i]);
  }

  TXTYK res = fitHelixSecondOrder(meas);
  printf("tx=%e\n",res.tx);
  printf("ty=%e\n",res.ty);  
  printf("k=%e\n",res.k);

}
