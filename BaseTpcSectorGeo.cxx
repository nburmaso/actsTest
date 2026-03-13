//-----------------------------------------------------------
// Description:
//      BaseTpcSectorGeo base class source file
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Authors:
//      Alexander Bychkov, Slavomir Hnatic
//      JINR, October, 2022
//-----------------------------------------------------------

#include "BaseTpcSectorGeo.h"

//__________________________________________________________________________

BaseTpcSectorGeo::BaseTpcSectorGeo() {}

//__________________________________________________________________________

BaseTpcSectorGeo::~BaseTpcSectorGeo() {}

//__________________________________________________________________________

int BaseTpcSectorGeo::SectorNumberFromGlobal(const TVector3 &globalXYZ)
{
   double phiGlobalNormalized = TMath::ATan2(globalXYZ.X(), globalXYZ.Y()) / sectorPhiRad;
   // convention: 105 deg is start of Sector 0, 75 deg is start of Sector 1, ..
   const double PHI_SHIFT_NORMALIZED = 3.5;
   int          iSector;
   if (phiGlobalNormalized > PHI_SHIFT_NORMALIZED)
      iSector = static_cast<int>(PHI_SHIFT_NORMALIZED - phiGlobalNormalized + sectorCountHalf);
   else
      iSector = static_cast<int>(PHI_SHIFT_NORMALIZED - phiGlobalNormalized);

   if (globalXYZ.Z() < 0.0) iSector += sectorCountHalf;

   return iSector;
}

//__________________________________________________________________________

TVector2 BaseTpcSectorGeo::PadRow2Local(double pad, double row)
{
   double x, y;
   int    rowInt = static_cast<int>(row);
   if (row > rowCount[inner]) {
      y = yPadAreaLength[inner] + (row - rowCount[inner]) * padHeight[outer];
      x = (pad - padCount[rowInt]) * padWidth[outer];
   } else {
      y = row * padHeight[inner];
      x = (pad - padCount[rowInt]) * padWidth[inner];
   }

   return {x, y};
}

//__________________________________________________________________________

TVector2 BaseTpcSectorGeo::PadRowCenter2Local(int padNumber, int rowNumber)
{
   return PadRow2Local(padNumber + 0.5, rowNumber + 0.5);
}

//__________________________________________________________________________

std::pair<double, double> BaseTpcSectorGeo::Local2PadRow(const TVector3 &localXYZ)
{
   double                          pad, row;
   const std::pair<double, double> ROW_OUTOFRANGE(0, -1.);
   const std::pair<double, double> PAD_OUTOFRANGE(1e4, 0);

   PadArea currentPadArea;
   if (localXYZ.Y() >= yPadAreaLocal[lowerEdge] && localXYZ.Y() <= yPadAreaLocal[midBoundary]) {
      currentPadArea = inner;
      row            = localXYZ.Y() / padHeight[inner];
   } else if (localXYZ.Y() > yPadAreaLocal[midBoundary] && localXYZ.Y() < yPadAreaLocal[upperEdge]) {
      currentPadArea = outer;
      row            = rowCount[inner] + (localXYZ.Y() - yPadAreaLength[inner]) / padHeight[outer];
   } else
      return ROW_OUTOFRANGE;

   pad = localXYZ.X() / padWidth[currentPadArea] + padCount[static_cast<int>(row)];

   if (pad < 0 || pad >= 2 * padCount[static_cast<int>(row)]) return PAD_OUTOFRANGE;

   return std::make_pair(pad, row);
}

//__________________________________________________________________________

TVector3 BaseTpcSectorGeo::Global2Local(const TVector3 &globalXYZ, int iSector)
{
   iSector = SectorNumberFromGlobal(globalXYZ);

   TVector3 localXYZ(globalXYZ);

   localXYZ.RotateZ(sectorPhiRad * iSector);

   if (localXYZ.Z() > 0) localXYZ.RotateY(TMath::Pi());

   localXYZ.SetY(localXYZ.Y() - yPadAreaLowerEdge);
   localXYZ.SetZ(localXYZ.Z() + driftLength);

   return localXYZ;
}

//__________________________________________________________________________

std::pair<TVector3, int> BaseTpcSectorGeo::Global2Local(const TVector3 &globalXYZ)
{
   int      iSector;
   TVector3 localXYZ = Global2Local(globalXYZ, iSector);

   return std::make_pair(localXYZ, iSector);
}

//__________________________________________________________________________

TVector3 BaseTpcSectorGeo::Local2Global(const TVector3 &localXYZ, int iSector)
{

   TVector3 globalXYZ(localXYZ);

   globalXYZ.SetZ(localXYZ.Z() - driftLength);
   globalXYZ.SetY(localXYZ.Y() + yPadAreaLowerEdge);

   if (iSector < sectorCountHalf) globalXYZ.RotateY(TMath::Pi());

   globalXYZ.RotateZ(-sectorPhiRad * iSector);

   return globalXYZ;
}

//__________________________________________________________________________

TVector3 BaseTpcSectorGeo::Local2Global(const std::pair<TVector3, int> &localPosition)
{
   return Local2Global(localPosition.first, localPosition.second);
}
