void kinematics(){
  double p4 = 0.9;
  double eta4 = 2.2;
  double theta4 = 2 * TMath::ATan(TMath::Exp(-eta4));
  double pt4 = p4*TMath::Sin(theta4);
  printf("%f\n", pt4);
  return;


  double eta = 1.5;
  TVector3 v;
  v.SetPtEtaPhi(0.2,eta,0);
  printf("%f\n", tan(v.Theta()));
  printf("%f\n", v.Mag());


  double tanTheta = 1300/3000.;
  //double tanTheta = 340/1630.;
  TVector3 v2;
  v2.SetMagThetaPhi(0.2,atan(tanTheta),0);
  printf("%f\n", v2.Eta());

  double p3 = 0.4;
  double eta3 = 1.6;
  double theta3 = 2 * TMath::ATan(TMath::Exp(-eta3));
  double pt3 = p3*TMath::Sin(theta3);
  printf("%f\n", pt3);
  printf("%f\n",1./TMath::Tan(theta3));

}

