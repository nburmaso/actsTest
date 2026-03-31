#pragma once

#include <cmath>
#include <numbers>
#include <vector>

#include "TObject.h"

using std::numbers::pi;

class MyFtdGeo
{
public:
  MyFtdGeo() = default;
  virtual ~MyFtdGeo(){}

  enum FtdLayerTypes {
    kPixel = 2,
    kStrawR = 4,
    kStrawU = 5,
    kStrawV = 6
  };

  double GetField() const { return fBz; }

  double GetFtdZMin() const { return fFtdZMin; }
  double GetFtdZMax() const { return fFtdZMax; }

  double GetZMinEndCap() const { return fEndCapZMin; }
  double GetZMaxEndCap() const { return fEndCapZMax; }

  double GetRMax() const { return fRMax; }
  double GetLayerRMin(int iLayer) const { return fLayerRMin[iLayer]; }
  double GetLayerRMid(int iLayer) const { return fLayerRMid[iLayer]; }
  double GetLayerRMax(int iLayer) const { return fLayerRMax[iLayer]; }

  const std::vector<double>& GetLayerPositions() const { return fLayerPositions; }
  int    GetLayerNumberOfTubes(int iLayer) const { return fLayerNumberOfTubes[iLayer]; }
  double GetLayerPosition(int iLayer) const { return fLayerPositions[iLayer]; }
  double GetLayerThickness(int iLayer) const { return fLayerThickness[iLayer]; }
  double GetLayerAngle(int iLayer) const { return fLayerAngle[iLayer]; }
  int    GetLayerType(int iLayer) const { return fLayerType[iLayer]; }
  double GetLayerStripWidth(int iLayer) const { return fLayerStripWidth[iLayer]; }
  int    GetNumberOfLayers() const { return fLayerPositions.size(); }
  int    GetNumberOfStations() const { return fNStations; }
  int    GetNLayersPerStation() const { return fLayersPerStation; }
  int    GetLayerStation(int iLayer) const{ return fLayerStations[iLayer]; }
  double GetLayerStereoAngle(int iLayer) const { return fLayerType[iLayer]==5 ? fTubeIncl : (fLayerType[iLayer]==6 ? -fTubeIncl : 0); }
  double GetLayerCenterR(int iLayer) const { return (fLayerRMin[iLayer]+fLayerRMax[iLayer])/2.; }
  double GetTubeCenterAngle(int iLayer, int iTube) const { return fLayerAngle[iLayer] + 2*pi/fLayerNumberOfTubes[iLayer]*(iTube + 0.5); }
  double GetTubeCenterX(int iLayer, int iTube) const { return GetLayerCenterR(iLayer)*cos(GetTubeCenterAngle(iLayer, iTube)); }
  double GetTubeCenterY(int iLayer, int iTube) const { return GetLayerCenterR(iLayer)*sin(GetTubeCenterAngle(iLayer, iTube)); }
  double GetTubeRotationAngle(int iLayer, int iTube) const { return GetTubeCenterAngle(iLayer, iTube) + GetLayerStereoAngle(iLayer); }
  double GetTubeRotationAngleDeg(int iLayer, int iTube) const { return 180./pi*GetTubeRotationAngle(iLayer, iTube); }
  std::pair<double, double> GetTubeLocCoordinates(double x, double y, int iLayer, int iTube, bool debug = 0) {
    double phi = GetTubeRotationAngle(iLayer, iTube);
    double xc = GetTubeCenterX(iLayer, iTube);
    double yc = GetTubeCenterY(iLayer, iTube);
    double cosp = cos(phi);
    double sinp = sin(phi);
    double loc1 =  (x - xc) * cosp + (y - yc) * sinp;
    double loc0 = -(x - xc) * sinp + (y - yc) * cosp;
    if (debug) printf("phi=%f xc=%f yc=%f cos(phi)=%f sin(phi)=%f loc1=%f loc0=%f\n", phi, xc, yc, cosp, sinp, loc1, loc0);
    if (debug) printf("(x-xc)*sinp=%f (y-yc)*cosp=%f\n", -(x - xc) * sinp, (y - yc) * cosp);
    return {loc0, loc1};
  }

  double GetRocThick() const { return fRocThick; }
  double GetRocRMin() const { return fRocRMin; }
  double GetRocRMax() const { return fRocRMax; }

  double GetFrameThick() const { return fFrameThick; }
  double GetFrameRMin1() const { return fFrameRMin1; }
  double GetFrameRMax1() const { return fFrameRMax1; }
  double GetFrameRMin2() const { return fFrameRMin2; }
  double GetFrameRMax2() const { return fFrameRMax2; }

  double GetFrameHalfYRadial() const { return fFrameHalfYRadial; }
  double GetFrameHalfXRadial() const { return fFrameHalfXRadial; }
  double GetFrameCXRadial0() const { return fFrameCXRadial0; }

  double GetFtdRMin() const { return fFtdRMin; }
  double GetFtdRMax() const { return fFtdRMax; }

  double GetPipeR() const { return fPipeR; }
  double GetPipeThick() const { return fPipeThick; }

  double GetTubeInclDeg() const { return fTubeInclDeg; }
  double GetTubeIncl() const { return fTubeIncl; }

  bool IsGeometrySimple() const { return fGeometrySimple; }
private:
  const double fBz{0.5};  // T

  // 5 pixel-like FTD stations + 1 pixel-like TOF
//  const bool fGeometrySimple = kTRUE;
//  const std::vector<double> fLayerPositions{{210, 232.5, 255, 277.5, 300., 350.}}; // cm
//  const std::vector<int> fLayerType{{2, 2, 2, 2, 2, 2}};
//  const std::vector<double> fLayerThickness{{0.01, 0.01, 0.01, 0.01, 0.01, 0.01}}; // cm
//  const std::vector<double> fLayerAngle{{0., 0., 0., 0., 0., 0.}}; // rad dummy
//  const std::vector<double> fLayerStripWidth{{0.5, 0.5, 0.5, 0.5, 0.5, 0.5}}; // cm dummy
//  const std::vector<double> fLayerRMin{{61, 67, 74, 80, 87, 101}}; // cm
//  const std::vector<double> fLayerRMid{{77, 85, 94, 102, 111, 126}}; // cm dummy
//  const std::vector<double> fLayerRMax{{94, 104, 114, 124, 134, 156}}; // cm
//  const std::vector<int>    fLayerNumberOfTubes{{1, 1, 1, 1, 1, 1}}; // dummy

  // 4 strip-like stations + 1 pixel-like station
  // const std::vector<double> fLayerPositions{{210, 211, 232, 233, 255, 256, 277, 278, 300}}; // cm
  // const std::vector<int> fLayerType{{1, 1, 1, 1, 1, 1, 1, 1, 2}};
  // const std::vector<double> fLayerThickness{{0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01}}; // cm
  // const std::vector<double> fLayerAngle     {{0., pi/2, 0., pi/2, 0., pi/2, 0., pi/2, 0.}}; // rad
  // const std::vector<double> fLayerStripWidth{{0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}}; // cm
  // const std::vector<double> fLayerRMin{{61, 61, 67, 67, 74, 74, 80, 80, 87}}; // cm
  // const std::vector<double> fLayerRMid{{77, 77, 85, 85, 94, 94, 102, 102, 111}}; // cm
  // const std::vector<double> fLayerRMax{{94, 94, 104, 104, 114, 114, 124, 124, 134}}; // cm  
  // const std::vector<int>    fLayerNumberOfTubes{{1, 1, 1, 1, 1, 1, 1, 1, 1}}; // dummy

  //4 strip-like stations + 1 pixel-like station 1 cm 1 ring
  // const std::vector<double> fLayerPositions{{210, 211, 232, 233, 255, 256, 277, 278, 300}}; // cm
  // const std::vector<int> fLayerType{{1, 1, 1, 1, 1, 1, 1, 1, 2}};
  // const std::vector<double> fLayerThickness{{0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01}}; // cm
  // const std::vector<double> fLayerAngle     {{0., pi/2, 0., pi/2, 0., pi/2, 0., pi/2, 0.}}; // rad
  // const std::vector<double> fLayerStripWidth{{1., 1., 1., 1., 1., 1., 1., 1., 1.}}; // cm
  // const std::vector<double> fLayerRMin{{60, 60, 60, 60, 60, 60, 60, 60, 60}}; // cm
  // const std::vector<double> fLayerRMid{{130, 130, 130, 130, 130, 130, 130, 130, 130}}; // cm
  // const std::vector<double> fLayerRMax{{131, 131, 131, 131, 131, 131, 131, 131, 131}}; // cm

  // 4 strip-like stations + 1 pixel-like station 0.5 cm 1 small ring
  // const std::vector<double> fLayerPositions{{210, 211, 232, 233, 255, 256, 277, 278, 300}}; // cm
  // const std::vector<int> fLayerType{{1, 1, 1, 1, 1, 1, 1, 1, 2}};
  // const std::vector<double> fLayerThickness{{0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01}}; // cm
  // const std::vector<double> fLayerAngle     {{0., pi/2, 0., pi/2, 0., pi/2, 0., pi/2, 0.}}; // rad
  // const std::vector<double> fLayerStripWidth{{0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}}; // cm
  // const std::vector<double> fLayerRMin{{61, 61, 67, 67, 74, 74, 80, 80, 87}}; // cm
  // const std::vector<double> fLayerRMid{{93, 93, 103, 103, 113, 113, 123, 123, 133}}; // cm
  // const std::vector<double> fLayerRMax{{94, 94, 104, 104, 114, 114, 124, 124, 134}}; // cm  

  // 5 pixel-like station
  // const std::vector<double> fLayerPositions{{210, 232.5, 255, 277.5, 300}}; // cm
  // const std::vector<int> fLayerType{{2, 2, 2, 2, 2}};
  // const std::vector<double> fLayerThickness{{0.01, 0.01, 0.01, 0.01, 0.01}}; // cm
  // const std::vector<double> fLayerAngle{{0., 0., 0., 0., 0.}}; // rad
  // const std::vector<double> fLayerStripWidth{{0.5, 0.5, 0.5, 0.5, 0.5}}; // cm
  // const std::vector<double> fLayerRMin{{61, 67, 74, 80, 87}}; // cm
  // const std::vector<double> fLayerRMid{{77, 85, 94, 102, 111}}; // cm
  // const std::vector<double> fLayerRMax{{94, 104, 114, 124, 134}}; // cm  

  // 5 stations x 6 layers, 4 mm RUV straw-tubes
  // const std::vector<double> fLayerPositions{{210, 211, 232, 233, 255, 256, 277, 278, 300}}; // cm
  // const std::vector<int> fLayerType{{1, 1, 1, 1, 1, 1, 1, 1, 1}};
  // const std::vector<double> fLayerThickness{{0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6}}; // cm
  // const std::vector<double> fLayerAngle     {{0., 0.0082, 0., 0.0082, 0., 0.0082, 0., 0.0082, 0.}}; // rad
  // const std::vector<double> fLayerStripWidth{{1., 1., 1., 1., 1., 1., 1., 1., 1.}}; // cm
  // const std::vector<double> fLayerRMin{{61, 61, 61, 61, 61, 61, 61, 61, 61}}; // cm
  // const std::vector<double> fLayerRMid{{130, 130, 130, 130, 130, 130, 130, 130, 130}}; // cm
  // const std::vector<double> fLayerRMax{{93, 93, 93, 93, 93, 93, 93, 93, 93}}; // cm
  // const std::vector<int>    fLayerNumberOfTubes{{766, 766, 766, 766, 766, 766, 766, 766, 766}};

  // const std::vector<double> fLayerPositions{{210, 211}}; // cm
  // const std::vector<int> fLayerType{{4, 4}};
  // const std::vector<double> fLayerThickness{{0.6, 0.6}}; // cm
  // const std::vector<double> fLayerAngle     {{0., 0.0041}}; // rad
  // const std::vector<double> fLayerStripWidth{{0.4, 0.4}}; // cm
  // const std::vector<double> fLayerRMin{{61, 61}}; // cm
  // const std::vector<double> fLayerRMid{{130, 130}}; // cm
  // const std::vector<double> fLayerRMax{{93, 93}}; // cm
  // const std::vector<int>    fLayerNumberOfTubes{{766, 766}};

  // const std::vector<double> fLayerPositions{{210, 210.8, 211.6}}; // cm
  // const std::vector<int> fLayerType{{4, 5, 6}};
  // const std::vector<double> fLayerThickness{{0.6, 0.6, 0.6}}; // cm
  // const std::vector<double> fLayerAngle     {{0., 0., 0.}}; // rad
  // const std::vector<double> fLayerStripWidth{{0.4, 0.4, 0.4}}; // cm
  // const std::vector<double> fLayerRMin{{61, 61, 61.}}; // cm
  // const std::vector<double> fLayerRMid{{130, 130, 130}}; // cm
  // const std::vector<double> fLayerRMax{{93, 93, 93}}; // cm
  // const std::vector<int>    fLayerNumberOfTubes{{766, 766, 766}};

  // // 5 stations with 6 layers of 4mm straw tubes + 1 fake pixel-like layer for seeding
  // const bool fGeometrySimple = false;
  // const int fLayersPerStation = 7;
  // const std::vector<double> fLayerPositions{{210, 210.8, 211.6, 212., 212.4, 213.2, 214.,
  //                                            232, 232.8, 233.6, 234., 234.4, 235.2, 236.,
  //                                            254, 254.8, 255.6, 256., 256.4, 257.2, 258.,
  //                                            276, 276.8, 277.6, 278., 278.4, 279.2, 280.,
  //                                            298, 298.8, 299.6, 300., 300.4, 301.2, 302.
  // }}; // cm

  // const std::vector<int> fLayerType{{4, 5, 6, 2, 4, 5, 6,
  //                                    4, 5, 6, 2, 4, 5, 6,
  //                                    4, 5, 6, 2, 4, 5, 6,
  //                                    4, 5, 6, 2, 4, 5, 6,
  //                                    4, 5, 6, 2, 4, 5, 6
  // }};

  // const std::vector<double> fLayerThickness{{0.004, 0.004, 0.004, 0.0001, 0.004, 0.004, 0.004,
  //                                            0.004, 0.004, 0.004, 0.0001, 0.004, 0.004, 0.004,
  //                                            0.004, 0.004, 0.004, 0.0001, 0.004, 0.004, 0.004,
  //                                            0.004, 0.004, 0.004, 0.0001, 0.004, 0.004, 0.004,
  //                                            0.004, 0.004, 0.004, 0.0001, 0.004, 0.004, 0.004
  // }}; // cm

  // const std::vector<double> fLayerAngle     {{0., 0., 0., 0., 0.0041, 0.0041, 0.0041,
  //                                            0., 0., 0., 0., 0.0041, 0.0041, 0.0041,
  //                                            0., 0., 0., 0., 0.0041, 0.0041, 0.0041,
  //                                            0., 0., 0., 0., 0.0041, 0.0041, 0.0041,
  //                                            0., 0., 0., 0., 0.0041, 0.0041, 0.0041
  // }}; // rad

  // const std::vector<double> fLayerStripWidth{{0.4, 0.4, 0.4, 0., 0.4, 0.4, 0.4,
  //                                            0.4, 0.4, 0.4, 0., 0.4, 0.4, 0.4,
  //                                            0.4, 0.4, 0.4, 0., 0.4, 0.4, 0.4,
  //                                            0.4, 0.4, 0.4, 0., 0.4, 0.4, 0.4,
  //                                            0.4, 0.4, 0.4, 0., 0.4, 0.4, 0.4
  // }}; // cm

  // const std::vector<double> fLayerRMin{{61., 61., 61., 61., 61., 61., 61.,
  //                                      67., 67., 67., 67., 67., 67., 67.,
  //                                      73., 73., 73., 73., 73., 73., 73.,
  //                                      80., 80., 80., 80., 80., 80., 80.,
  //                                      86., 86., 86., 86., 86., 86., 86.
  // }}; // cm

  // const std::vector<double> fLayerRMax{{  93., 93., 93., 93., 93., 93., 93.,
  //                                        103.,103.,103.,103.,103.,103.,103.,
  //                                        113.,113.,113.,113.,113.,113.,113.,
  //                                        123.,123.,123.,123.,123.,123.,123.,
  //                                        133.,133.,133.,133.,133.,133.,133.
  // }}; // cm

  // const std::vector<int>fLayerNumberOfTubes{{766, 766, 766, 0, 766, 766, 766,
  //                                            846, 846, 846, 0, 846, 846, 846,
  //                                            927, 927, 927, 0, 927, 927, 927,
  //                                           1007,1007,1007, 0,1007,1007,1007,
  //                                           1087,1087,1087, 0,1087,1087,1087
  // }};

  // const std::vector<double> fLayerRMid{{130, 130, 130., 130, 130, 130, 130.,
  //                                       130, 130, 130., 130, 130, 130, 130.,
  //                                       130, 130, 130., 130, 130, 130, 130.,
  //                                       130, 130, 130., 130, 130, 130, 130.,
  //                                       130, 130, 130., 130, 130, 130, 130.,
  // }}; // cm

// 5 stations with 6 layers of 5mm straw tubes + 1 fake pixel-like layer for seeding
  // const bool fGeometrySimple = false;
  // const int fNStations = 5;
  // const int fLayersPerStation = 7;

  // const std::vector<int> fLayerStations{{
  //   0,
  //   0, 0, 0, 0, 0, 0, 0,
  //   1, 1, 1, 1, 1, 1, 1,
  //   2, 2, 2, 2, 2, 2, 2,
  //   3, 3, 3, 3, 3, 3, 3,
  //   4, 4, 4, 4, 4, 4, 4
  //  ,4
  // }}; // cm

  // const std::vector<double> fLayerPositions{{
  //   207,
  //   210, 210.8, 211.6, 212., 212.4, 213.2, 214.,
  //   232, 232.8, 233.6, 234., 234.4, 235.2, 236.,
  //   254, 254.8, 255.6, 256., 256.4, 257.2, 258.,
  //   276, 276.8, 277.6, 278., 278.4, 279.2, 280.,
  //   298, 298.8, 299.6, 300., 300.4, 301.2, 302.
  //  ,305
  // }}; // cm

  // const std::vector<int> fLayerType{{
  //   2,
  //   4, 5, 6, 2, 4, 5, 6,
  //   4, 5, 6, 2, 4, 5, 6,
  //   4, 5, 6, 2, 4, 5, 6,
  //   4, 5, 6, 2, 4, 5, 6,
  //   4, 5, 6, 2, 4, 5, 6, 2
  // }};

  // const std::vector<double> fLayerThickness{{
  //   0.0001,
  //   0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025,
  //   0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025,
  //   0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025,
  //   0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025,
  //   0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025, 0.0001
  // }}; // cm

  // const std::vector<double> fLayerAngle     {{
  //   0.,
  //   0., 0., 0., 0., 0.0052, 0.0052, 0.0052,
  //   0., 0., 0., 0., 0.0047, 0.0047, 0.0047,
  //   0., 0., 0., 0., 0.0043, 0.0043, 0.0043,
  //   0., 0., 0., 0., 0.0039, 0.0039, 0.0039,
  //   0., 0., 0., 0., 0.0036, 0.0036, 0.0036, 0.
  // }}; // rad

  // const std::vector<double> fLayerStripWidth{{
  //   0.,
  //   0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5,
  //   0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5,
  //   0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5,
  //   0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5,
  //   0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5, 0.
  // }}; // cm

  // const std::vector<double> fLayerRMin{{
  //   58.,
  //   58.,  58.,  58.,  58.,  58.,  58., 58.,
  //   64.,  64.,  64.,  64.,  64.,  64., 64.,
  //   70.,  70.,  70.,  70.,  70.,  70., 70.,
  //   76.,  76.,  76.,  76.,  76.,  76., 76.,
  //   82.,  82.,  82.,  82.,  82.,  82., 82., 82.
  // }}; // cm

  // const std::vector<double> fLayerRMax{{
  //    95.,
  //    95.,  95.,  95.,  95.,  95.,  95.,  95.,
  //   105., 105., 105., 105., 105., 105., 105.,
  //   115., 115., 115., 115., 115., 115., 115.,
  //   124., 124., 124., 124., 124., 124., 124.,
  //   134., 134., 134., 134., 134., 134., 134., 134.
  // }}; // cm

  // const std::vector<int> fLayerNumberOfTubes{{
  //      0,
  //    606,  606,  606,    0,  606,  606, 606,
  //    669,  669,  669,    0,  669,  669, 669,
  //    733,  733,  733,    0,  733,  733, 733,
  //    796,  796,  796,    0,  796,  796, 796,
  //    860,  860,  860,    0,  860,  860, 860, 0
  // }};

  // const std::vector<double> fLayerRMid{{
  //   130,
  //   130, 130, 130., 130, 130, 130, 130.,
  //   130, 130, 130., 130, 130, 130, 130.,
  //   130, 130, 130., 130, 130, 130, 130.,
  //   130, 130, 130., 130, 130, 130, 130.,
  //   130, 130, 130., 130, 130, 130, 130., 130.
  // }}; // cm

// 5 stations with 8 layers of 5mm straw tubes + 1 fake pixel-like layer
  const bool fGeometrySimple = false;
  const int fNStations = 5;
  const int fLayersPerStation = 9;

  const std::vector<int> fLayerStations{{
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 1, 1, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 2, 2, 2, 2,
    3, 3, 3, 3, 3, 3, 3, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4, 4
  }}; // cm

  const std::vector<double> fLayerPositions{{
    210, 210.7, 211.4, 212.1, 212.4, 212.8, 213.5, 214.2, 214.9,
    232, 232.7, 233.4, 234.1, 234.4, 234.8, 235.5, 236.2, 236.9,
    254, 254.7, 255.4, 256.1, 256.4, 256.8, 257.5, 258.2, 258.9,
    276, 276.7, 277.4, 278.1, 278.4, 278.8, 279.5, 280.2, 280.9,
    298, 298.7, 299.4, 300.1, 300.4, 300.8, 301.5, 302.2, 302.9
  }}; // cm

  const std::vector<int> fLayerType{{
    5, 5, 6, 6, 2, 5, 5, 6, 6,
    5, 5, 6, 6, 2, 5, 5, 6, 6,
    5, 5, 6, 6, 2, 5, 5, 6, 6,
    5, 5, 6, 6, 2, 5, 5, 6, 6,
    5, 5, 6, 6, 2, 5, 5, 6, 6
  }};

  const std::vector<double> fLayerThickness{{
    0.0025, 0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025, 0.0025,
    0.0025, 0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025, 0.0025,
    0.0025, 0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025, 0.0025,
    0.0025, 0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025, 0.0025,
    0.0025, 0.0025, 0.0025, 0.0025, 0.0001, 0.0025, 0.0025, 0.0025, 0.0025
  }}; // cm

  const std::vector<double> fLayerAngle {{
    0., 0.00260, 0., 0.00260, 0., 0.00520, 0.00780, 0.00520, 0.00780,
    0., 0.00235, 0., 0.00235, 0., 0.00470, 0.00705, 0.00470, 0.00705,
    0., 0.00215, 0., 0.00215, 0., 0.00430, 0.00645, 0.00430, 0.00645,
    0., 0.00195, 0., 0.00195, 0., 0.00390, 0.00585, 0.00390, 0.00585,
    0., 0.00183, 0., 0.00183, 0., 0.00365, 0.00548, 0.00365, 0.00548
  }}; // rad

  const std::vector<double> fLayerStripWidth{{
    0.5, 0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5, 0.5,
    0.5, 0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5, 0.5,
    0.5, 0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5, 0.5,
    0.5, 0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5, 0.5,
    0.5, 0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5, 0.5
  }}; // cm

  const std::vector<double> fLayerRMin{{
    58.,  58.,  58.,  58.,  58.,  58., 58., 58., 58.,
    64.,  64.,  64.,  64.,  64.,  64., 64., 64., 64.,
    70.,  70.,  70.,  70.,  70.,  70., 70., 70., 70.,
    76.,  76.,  76.,  76.,  76.,  76., 76., 76., 76.,
    82.,  82.,  82.,  82.,  82.,  82., 82., 82., 82.
  }}; // cm

  const std::vector<double> fLayerRMax{{
     95.,  95.,  95.,  95.,  95.,  95.,  95.,  95.,  95.,
    105., 105., 105., 105., 105., 105., 105., 105., 105.,
    115., 115., 115., 115., 115., 115., 115., 115., 115.,
    124., 124., 124., 124., 124., 124., 124., 124., 124.,
    134., 134., 134., 134., 134., 134., 134., 134., 134.,
  }}; // cm

  const std::vector<int> fLayerNumberOfTubes{{
     606, 606,  606,  606,    0,  606,  606, 606, 606,
     669, 669,  669,  669,    0,  669,  669, 669, 669,
     733, 733,  733,  733,    0,  733,  733, 733, 733,
     796, 796,  796,  796,    0,  796,  796, 796, 796,
     860, 860,  860,  860,    0,  860,  860, 860, 860
  }};

  const std::vector<double> fLayerRMid{{
    130, 130, 130., 130, 130, 130, 130., 130, 130.,
    130, 130, 130., 130, 130, 130, 130., 130, 130.,
    130, 130, 130., 130, 130, 130, 130., 130, 130.,
    130, 130, 130., 130, 130, 130, 130., 130, 130.,
    130, 130, 130., 130, 130, 130, 130., 130, 130.
  }}; // cm

  const double fEps = 1e-12;

  // straw tube parameters
  const double fTubeInclDeg{7.};
  const double fTubeIncl{fTubeInclDeg / 180. * pi};

  // toy pipe parameters
  const double fPipeRMin{3.1}; // cm
  const double fPipeRMax{3.225}; // cm
  const double fPipeThick{fPipeRMax - fPipeRMin}; // cm
  const double fPipeR{0.5 * (fPipeRMax + fPipeRMin)}; // cm

  // ftd volume boundaries
  const double fFtdZMin{fLayerPositions.front() - 0.5}; // cm
  const double fFtdZMax{fLayerPositions.back() + 0.5}; // cm

  // roc parameters
  const double fRocThick{0.25 * 8.897}; // cm
  const double fRocRMin{27.}; // cm
  const double fRocRMax{141.}; // cm

  // frame parameters
  const double fFrameThick{0.85 * 8.897}; // cm
  const double fFrameRMin1{35.}; // cm
  const double fFrameRMax1{42.}; // cm
  const double fFrameRMin2{120.}; // cm
  const double fFrameRMax2{141.}; // cm
  const double fFrameHalfYRadial{4.}; // cm
  const double fFrameHalfXRadial{(std::sqrt(fFrameRMin2 * fFrameRMin2 - fFrameHalfYRadial * fFrameHalfYRadial) - fFrameRMax1) / 2. - 2. * fEps}; // cm
  const double fFrameCXRadial0{fFrameRMax1 + fFrameHalfXRadial + fEps}; // cm

  // endcap volume boundaries
  const double fEndCapZMin{165.}; // cm // todo: 163.879 + 0.01 + [pad thickness?]
  const double fEndCapZMax{fLayerPositions.front() - 0.2}; // cm
  const double fRMax{std::sqrt(2.) * fLayerRMax.back() + 0.2}; // ~ 184 cm
  const double fFtdRMin{fFrameRMin1 - fEps}; // cm
  const double fFtdRMax{fFrameRMax2 + fEps}; // cm
};
