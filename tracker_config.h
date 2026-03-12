#pragma once

double bz = 0.5;    // T

//std::vector<double> positions{{210, 232.5, 255, 277.5, 300}}; // z-positions of sensitive layers in cm
//std::vector<double> positions{{210, 211, 232.5, 233.5, 255, 256, 277.5, 278.5, 300}}; // double layers for strip-like simulations
const std::vector<double> positions{{
    207,
    210, 210.8, 211.6, 212., 212.4, 213.2, 214.,
    232, 232.8, 233.6, 234., 234.4, 235.2, 236.,
    254, 254.8, 255.6, 256., 256.4, 257.2, 258.,
    276, 276.8, 277.6, 278., 278.4, 279.2, 280.,
    298, 298.8, 299.6, 300., 300.4, 301.2, 302.
   ,305
}}; // cm

const std::vector<double> layerAngle {{
    0.,
    0., 0., 0., 0., 0.0052, 0.0052, 0.0052,
    0., 0., 0., 0., 0.0047, 0.0047, 0.0047,
    0., 0., 0., 0., 0.0043, 0.0043, 0.0043,
    0., 0., 0., 0., 0.0039, 0.0039, 0.0039,
    0., 0., 0., 0., 0.0036, 0.0036, 0.0036
   ,0.
}}; // rad

const std::vector<double> layerRMin {{
     58.,
     58.,  58.,  58.,  58.,  58.,  58., 58.,
     64.,  64.,  64.,  64.,  64.,  64., 64.,
     70.,  70.,  70.,  70.,  70.,  70., 70.,
     76.,  76.,  76.,  76.,  76.,  76., 76.,
     82.,  82.,  82.,  82.,  82.,  82., 82.
    ,82.
}}; // cm

const std::vector<double> layerRMax{{
    95.,
    95.,  95.,  95.,  95.,  95.,  95.,  95.,
    105., 105., 105., 105., 105., 105., 105.,
    115., 115., 115., 115., 115., 115., 115.,
    124., 124., 124., 124., 124., 124., 124.,
    134., 134., 134., 134., 134., 134., 134.
   ,134.
}}; // cm

const std::vector<int> numberOfTubes{{
       0,
     606,  606,  606,    0,  606,  606, 606,
     669,  669,  669,    0,  669,  669, 669,
     733,  733,  733,    0,  733,  733, 733,
     796,  796,  796,    0,  796,  796, 796,
     860,  860,  860,    0,  860,  860, 860
      ,0
}};

const std::vector<int> layerType{{
    2,
    4, 5, 6, 2, 4, 5, 6,
    4, 5, 6, 2, 4, 5, 6,
    4, 5, 6, 2, 4, 5, 6,
    4, 5, 6, 2, 4, 5, 6,
    4, 5, 6, 2, 4, 5, 6
   ,2
}};

const std::vector<double> layerThickness{{
    0.0001,
    0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025,
    0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025,
    0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025,
    0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025,
    0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025
   ,0.0001
}}; // cm


double angleRotDeg = 7.; // degrees
double layerStripWidth = 0.5; // cm

double radLenSilicon = 9.370;  // cm
double radLenAluminium = 8.897; // cm

double rMinStation = 35.7; // cm
double rMaxStation = 135.; // cm

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
