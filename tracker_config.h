#pragma once

double bz = 0.5;    // T

std::vector<double> positions{{210, 232.5, 255, 277.5, 300}}; // z-positions of sensitive layers in cm
double radLenSilicon = 9.370;  // cm
double radLenAluminium = 8.897; // cm

// double thickness{0.00112*radLenSilicon};
double thickness = 0.01; // cm
double rMinStation = 35.7; // cm
double rMaxStation = 130.; // cm

double eps = 1e-10; // cm
double radLenFractionROC = 0.25;
double thicknessROC = radLenFractionROC*radLenAluminium;
double zDrift = 163; // cm
double zROC = zDrift + thicknessROC/2.;
double rMinROC = 27; // cm
double rMaxROC = 141; // cm
 
double radLenFractionFrame = 0.85;
double thicknessFrame = radLenFractionFrame*radLenAluminium;

double zFrame = zDrift + thicknessROC + thicknessFrame/2. + 0.1; // cm
double zFrameCircum1 = zFrame;
double rMinFrameCircum1 = 35; // cm
double rMaxFrameCircum1 = 42; // cm

double zFrameCircum2 = zFrame;
double rMinFrameCircum2 = 120; // cm
double rMaxFrameCircum2 = 141; // cm

double zFrameRadial = zFrame;
double halfYFrameRadial = 4.0;  // cm
double halfXFrameRadial = (sqrt(rMinFrameCircum2*rMinFrameCircum2-halfYFrameRadial*halfYFrameRadial)-rMaxFrameCircum1)/2.-2*eps;
double cXFrameRadial0 = rMaxFrameCircum1 + halfXFrameRadial + eps;
const int nSectors = 12;


double fwdRMin = 0.;
double fwdRMax = rMaxFrameCircum2 + eps;
double fwdHalfZ = positions.back()+10; // cm
