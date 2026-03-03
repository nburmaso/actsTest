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
    const int n = z.size();
    if (c.size() != n || s.size() != n || g.size() != n)
        throw std::runtime_error("Input vectors must have same length");

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

#include <iostream>

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
  double sigma = 0.001;

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
  printf("k=%f\n",res.k);

}
