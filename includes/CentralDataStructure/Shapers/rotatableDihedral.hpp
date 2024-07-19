#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ROTATABLEDIHEDRAL_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ROTATABLEDIHEDRAL_HPP

#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngles.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

// This class stores the four atoms that define a dihedral angle, the atoms that move when it is rotated
// and, if moved, the previous dihedral angle, which allows me to reset easily.
namespace cds
{
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

    struct RotatableDihedral
    {
        RotatableDihedral(bool isBranchingLinkage, const std::array<Atom*, 4>& atoms_,
                          const std::vector<DihedralAngleData>& metadata)
            : isBranchingLinkage(isBranchingLinkage), atoms(atoms_), metadataVector(metadata) {};

        bool isBranchingLinkage;
        // The four atoms that define the dihedral angle. The bond between atom2_ and atom3_ is what is rotated.
        std::array<cds::Atom*, 4> atoms;
        // A vector of pointers to the atoms that are connected to atom2_ and atom3_, and will be rotated when that bond
        // is rotated.
        std::vector<cds::Coordinate*> movingCoordinates;
        DihedralAngleDataVector metadataVector;
        DihedralAngleData currentMetadata;
    };

    std::string print(const RotatableDihedral& dihedral);
    DihedralCoordinates dihedralCoordinates(const cds::RotatableDihedral& dihedral);
    void setDihedralAngle(RotatableDihedral& dihedral, AngleWithMetadata target);
    bool setSpecificShape(RotatableDihedral& dihedral, std::string dihedralName, std::string selectedRotamer);

} // namespace cds
#endif
