#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_CDSGEOMETRYTOPOLOGY_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_CDSGEOMETRYTOPOLOGY_HPP_

#include "includes/CentralDataStructure/coordinate.hpp"

// Everything here should be ignorant of cds classes higher than coordinate.
namespace cds
{
double** GenerateRotationMatrix(Coordinate* direction, Coordinate* parent, double angle);
void SetDihedralAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, Coordinate* a4, const double dihedral_angle, std::vector<Coordinate*>& movingCoords);
void SetAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, const double angle, std::vector<Coordinate*> coordinatesToMove);

} // namespace

#endif /* INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_CDSGEOMETRYTOPOLOGY_HPP_ */