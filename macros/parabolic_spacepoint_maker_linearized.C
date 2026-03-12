#include <Eigen/Dense>
#include <vector>
#include <stdexcept>

struct HelixFitResult {
    double tx;
    double ty;
    double k;
};

HelixFitResult fitHelixLinear(
    const std::vector<double>& z,
    const std::vector<double>& c,
    const std::vector<double>& s,
    const std::vector<double>& g
) {
    const int n = c.size();
    // if (c.size() != n || s.size() != n || g.size() != n)
    //     throw std::runtime_error("Input vectors must have same length");

    Eigen::MatrixXd A(n, 4);
    Eigen::VectorXd b(n);

    for (int i = 0; i < n; ++i) {
        A(i, 0) =  s[i];           // t_x
        A(i, 1) = -c[i];           // t_y
        A(i, 2) =  z[i] * c[i];    // u_x
        A(i, 3) =  z[i] * s[i];    // u_y

        b(i) = -g[i];
    }

    // Solve least squares: min ||A*theta - b||
    Eigen::VectorXd theta =
        (A.transpose() * A).ldlt().solve(A.transpose() * b);

    double tx = theta(0);
    double ty = theta(1);
    double ux = theta(2);
    double uy = theta(3);

    double denom = tx * tx + ty * ty;
    if (denom == 0.0)
        throw std::runtime_error("Cannot extract curvature (zero slope)");

    double k = (ux * tx + uy * ty) / denom;

    return {tx, ty, k};
}

void parabolic_spacepoint_maker_linearized() {
  const int n = 5;
  std::vector<double> r = { 100,  100,  100,  100,  100}; // station minimal radius
  std::vector<double> z = {2100, 2110, 2120, 2130, 2140}; // station pozition
  double thetaTrueDeg = 30;
  double alphaTrueDeg = 20;
  double rTrue = 3000;
  double thetaTrue = thetaTrueDeg*TMath::DegToRad();
  double alphaTrue = alphaTrueDeg*TMath::DegToRad();
  std::vector<double> phiDeg = {alphaTrueDeg-4., alphaTrueDeg+1., alphaTrueDeg+6., alphaTrueDeg-3., alphaTrueDeg+3.};

  std::vector<double> d(n);
  std::vector<double> c(n);
  std::vector<double> s(n);
  std::vector<double> g(n);
  double sigma = 0.01;

  double tTrue = tan(thetaTrue);
  double pTrue = tTrue/rTrue; 
  std::vector<double> phi(n);
  std::vector<double> xTrue(n);
  std::vector<double> yTrue(n);
  std::vector<double> xmeas(n);
  std::vector<double> ymeas(n);
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
  auto res = fitHelixLinear(z,c,s,g);
  printf("tx=%f\n",res.tx);
  printf("ty=%f\n",res.ty);
  printf("k=%e\n",res.k);

  double chi2 = 0;
  for (int i=0;i<n;i++){
    double chi = (s[i]*res.tx - c[i]*res.ty + z[i]*c[i]*res.k*res.tx + z[i]*s[i]*res.k*res.ty + g[i])/d[i];
    chi2 += chi*chi;
  }
  printf("chi2=%f\n",chi2);


  // minimization algorithm  
  double xx[3];
  xx[0] = kTrue;
  xx[1] = txTrue;
  xx[2] = tyTrue;

  auto func = [&n,&s,&c,&z,&g](const double* xx) {
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
  };

  printf("%f\n",func(xx));

  ROOT::Math::Functor f(func,3);
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
//   kk = xs[0];
//   tx = xs[1];
//   ty = xs[2];
//   printf("tx=%f\n",tx);
//   printf("ty=%f\n",ty);
//   printf("k=%e\n",kk);
  double chi2min = 0;
  for (int i=0;i<n;i++){
    double chi = (s[i]*xs[1] - c[i]*xs[2] + z[i]*c[i]*xs[0]*xs[1] + z[i]*s[i]*xs[0]*xs[2] + g[i])/d[i];
    chi2min += chi*chi;
  }
  printf("chi2min=%f\n",chi2min);

}
