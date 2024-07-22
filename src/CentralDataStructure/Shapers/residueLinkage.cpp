#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"

#include <functional>
#include <sstream>

using cds::ResidueLinkage;
using cds::RotatableDihedral;

namespace
{
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

    int
    getNumberOfShapesBase(std::function<const DihedralAngleDataVector(const DihedralAngleDataVector&)> selectedMetadata,
                          const std::vector<RotatableDihedral>&
                              dihedrals) // Can have conformers (sets of rotamers) or permutations of rotamers
    {
        if (dihedrals.empty() || dihedrals[0].metadataVector.empty())
        {
            return 1;
        }
        else
        {
            auto& dihedral   = dihedrals[0];
            auto rotamerType = dihedral.metadataVector[0].rotamer_type_;
            if (rotamerType == gmml::MolecularMetadata::GLYCAM::RotamerType::permutation)
            {
                int numberOfShapes = 1;
                for (auto& entry : dihedrals)
                {
                    numberOfShapes *= selectedMetadata(entry.metadataVector).size();
                }
                return numberOfShapes;
            }
            else if (rotamerType == gmml::MolecularMetadata::GLYCAM::RotamerType::conformer)
            { // Conformer should mean that each dihedral will have the same number of metadata entries.
                // numberOfShapes = RotatableDihedrals_.size(); // This was correct for ASN for the wrong reason. 4
                // conformers and 4 dihedrals...
                return selectedMetadata(dihedral.metadataVector).size();
            }
            else
            {
                std::string str =
                    "Error: Unknown rotamer type: " +
                    std::to_string(
                        static_cast<std::underlying_type_t<gmml::MolecularMetadata::GLYCAM::RotamerType>>(rotamerType));
                gmml::log(__LINE__, __FILE__, gmml::ERR, str);
                throw std::runtime_error(str);
            }
        }
    }
} // namespace

std::vector<RotatableDihedral>
cds::rotatableDihedralsWithMultipleRotamers(const std::vector<RotatableDihedral>& dihedrals)
{
    std::vector<RotatableDihedral> returningDihedrals;
    for (auto& entry : dihedrals)
    {
        if (entry.metadataVector.size() > 1)
        {
            returningDihedrals.push_back(entry);
        }
    }
    return returningDihedrals;
}

int cds::numberOfShapes(const std::vector<RotatableDihedral>& dihedrals)
{
    auto selectMetadata = [](const DihedralAngleDataVector& data)
    {
        return data;
    };
    return getNumberOfShapesBase(selectMetadata, dihedrals);
}

int cds::numberOfLikelyShapes(const std::vector<RotatableDihedral>& dihedrals)
{
    return getNumberOfShapesBase(likelyMetadata, dihedrals);
}

cds::DihedralCoordinates cds::dihedralCoordinates(const cds::RotatableDihedral& dihedral)
{
    auto& atoms                       = dihedral.atoms;
    std::array<Coordinate*, 4> coords = {atoms[3]->getCoordinate(), atoms[2]->getCoordinate(),
                                         atoms[1]->getCoordinate(), atoms[0]->getCoordinate()};
    return DihedralCoordinates {*coords[0], *coords[1], *coords[2], *coords[3]};
}

std::string cds::print(const ResidueLinkage& linkage)
{
    std::stringstream ss;
    auto& linkageResidues = linkage.link.residues;
    ss << "ResidueLinkage Index: " << linkage.index << ", Name: " << linkage.name
       << ", NumberOfShapes: " << numberOfShapes(linkage.rotatableDihedrals)
       << ", ids: " << linkageResidues.first->getStringId() << "@" << linkage.link.atoms.first->getName() << " -- "
       << linkageResidues.second->getStringId() << "@" << linkage.link.atoms.second->getName() << "\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    for (auto& rotatableDihedral : linkage.rotatableDihedrals)
    {
        ss << print(rotatableDihedral);
    }
    return ss.str();
}

std::string cds::print(const RotatableDihedral& dihedral)
{
    auto& atoms = dihedral.atoms;
    std::stringstream ss;
    ss << atoms[0]->getName() << ", " << atoms[1]->getName() << ", " << atoms[2]->getName() << ", "
       << atoms[3]->getName() << ": " << cds::angle(dihedralCoordinates(dihedral)) << ".\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    return ss.str();
}
