//-------------------------------------------------------------------------------------------------
// Description:
//      BaseTpcSectorGeo class is the base class for TPC Geometry,
//      with its configuration in ideal, non-misaligned state.
//      Descriptive drawing:
//      https://git.jinr.ru/nica/docs/-/blob/main/docs/mpdroot/coding/geometry/TPC_coordinates.svg
//
//      Any custom TPC geometry must inherit from this class.
//
//      Note: If you are overriding any method, add in front of it the 'virtual' keyword.
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Authors:
//      Alexander Bychkov, Slavomir Hnatic
//      JINR, October, 2022
//-------------------------------------------------------------------------------------------------

#ifndef BASETPCSECTORGEO_HH
#define BASETPCSECTORGEO_HH

// ROOT Class Headers ---------------
#include <TObject.h>
#include <TMath.h>
#include <TVector2.h>
#include <TVector3.h>

class BaseTpcSectorGeo : public TObject {
public:
   // Constructors/Destructors ---------
   BaseTpcSectorGeo();
   virtual ~BaseTpcSectorGeo();

   /* constants */
   enum PadArea : int { inner, outer };
   enum PadAreaBoundary : int { lowerEdge, midBoundary, upperEdge };

   /* getters */
   // NOTE: setters are FORBIDDEN in this class !!!
   int    GetSectorCount() const { return sectorCount; }         // total number of sectors
   int    GetSectorCountHalf() const { return sectorCountHalf; } // number of sectors on one side
   double GetSectorPhiRad() const { return sectorPhiRad; }       // sector angle
   double GetSectorPhi0Rad() const { return sectorPhi0Rad; }     // first sector angle

   double GetZMin() const { return zMin; }               // half of membrane thickness
   double GetZMax() const { return zMax; }               // distance to pads in cm
   double GetDriftLength() const { return driftLength; } // drift length

   int    GetTimeBinCount() const { return timeBinCount; }   // number of timebins
   double GetTimeBinLength() const { return timeBinLength; } // timebin length, ns
   double GetTimeDriftMax() const { return timeDriftMax; }

   const std::vector<int>    &GetRowCount() const { return rowCount; }   // number of padrows for inner & outer regions
   const std::vector<double> &GetPadHeight() const { return padHeight; } // height of pads in inner & outer regions
   const std::vector<double> &GetPadWidth() const { return padWidth; }   // width of pads in inner & outer regions
   const std::vector<int>    &GetPadCount() const { return padCount; }   // half number of pads in rows

   double GetYPadPlaneOffset() const { return yPadPlaneOffset; } // Y of padplane's low edge in global coordinates
   double GetYPadPlane2PadAreaOffset() const { return yPadPlane2PadAreaOffset; } // padplane's lower edge to padarea
   double GetYPadAreaLowerEdge() const { return yPadAreaLowerEdge; }

   // Y length of padareas: inner, outer
   const std::vector<double> &GetYPadAreaLength() const { return yPadAreaLength; }
   // Y coordinates of padareas in local & global coordinates: inner edge, boundary inner/outer, outer edge
   const std::vector<double> &GetYPadAreaLocal() const { return yPadAreaLocal; }
   const std::vector<double> &GetYPadAreaGlobal() const { return yPadAreaGlobal; }

   /* transformations */
   // sectors: angle and number
   double SectorAxisAngleRad(int iSector) { return iSector * sectorPhiRad; } // angle of sector axis in radians
   int    SectorNumberFromGlobal(const TVector3 &globalXYZ);

   // padplane: padrows to local coordinates and vice versa
   TVector2                  PadRow2Local(double pad, double row);
   TVector2                  PadRowCenter2Local(int padNumber, int rowNumber);
   std::pair<double, double> Local2PadRow(const TVector3 &localXYZ);

   // time and z axis conversions
   int    Time2TimeBin(double time) { return TMath::FloorNint(time / timeBinLength); }
   double Time2TimeBinDouble(double time) { return time / timeBinLength; }
   double TimeBin2Time(int timeBin) { return timeBin * timeBinLength; }
   double TimeBin2Z(int timeBin, double velocity) { return velocity * TimeBin2Time(timeBin); }
   int    Z2TimeBin(double z, double velocity) { return Time2TimeBin(z / velocity); }

   // global to local coordinates and vice versa
   TVector3                 Global2Local(const TVector3 &globalXYZ, int iSector);
   std::pair<TVector3, int> Global2Local(const TVector3 &globalXYZ);
   TVector3                 Local2Global(const TVector3 &localXYZ, int iSector);
   TVector3                 Local2Global(const std::pair<TVector3, int> &localPosition);

protected:
private:
   int    sectorCount     = 24;                     // total number of sectors
   int    sectorCountHalf = 12;                     // number of sectors on one side
   double sectorPhiRad    = 30 * TMath::DegToRad(); // sector angle
   double sectorPhi0Rad   = -sectorPhiRad / 2.;     // first sector angle

   double zMin        = 0.01;               // half of membrane thickness
   double driftLength = 163.879;            // drift length
   double zMax        = zMin + driftLength; // distance to pads in cm

   int    timeBinCount  = 310;  // number of timebins
   double timeBinLength = 100.; // timebin length, ns
   double timeDriftMax  = timeBinCount * timeBinLength;

   std::vector<int>    rowCount{27, 26};    // number of padrows for inner & outer regions
   std::vector<double> padHeight{1.2, 1.8}; // height of pads in inner & outer regions
   std::vector<double> padWidth{0.5, 0.5};  // width of pads in inner & outer regions
   std::vector<int>    padCount{20, 21, 21, 22, 23, 23, 24, 24, 25, 26, 26, 27, 28, 28, 29, 30, 30, 31,
                             32, 32, 33, 33, 34, 35, 35, 36, 37, 38, 39, 40, 41, 41, 42, 43, 44, 45,
                             46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62};
   // half number of pads in rows

   double yPadPlaneOffset         = 40.3; // Y of padplane's low edge in global coordinates
   double yPadPlane2PadAreaOffset = 0.4;  // distance between padplane's lower edge and padarea
   double yPadAreaLowerEdge       = yPadPlaneOffset + yPadPlane2PadAreaOffset;
   // Y length of padareas: inner, outer
   std::vector<double> yPadAreaLength{rowCount[inner] * padHeight[inner], rowCount[outer] * padHeight[outer]};
   // Y coordinates of padareas in local & global coordinates: inner edge, boundary inner/outer, outer edge
   std::vector<double> yPadAreaLocal{0., yPadAreaLength[inner], yPadAreaLength[inner] + yPadAreaLength[outer]};
   std::vector<double> yPadAreaGlobal{yPadAreaLowerEdge, yPadAreaLowerEdge + yPadAreaLocal[midBoundary],
                                      yPadAreaLowerEdge + yPadAreaLocal[upperEdge]};

//   ClassDef(BaseTpcSectorGeo, 1);
};

#endif
