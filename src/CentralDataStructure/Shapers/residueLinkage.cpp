#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <functional>
#include <sstream>

using cds::ResidueLinkage;
using cds::RotatableDihedral;

namespace
{
    using GlycamMetadata::DihedralAngleDataVector;

    int
    getNumberOfShapesBase(std::function<const DihedralAngleDataVector(const DihedralAngleDataVector&)> selectedMetadata,
                          GlycamMetadata::RotamerType rotamerType, const std::vector<DihedralAngleDataVector>& metadata)
    {
        if (rotamerType == GlycamMetadata::RotamerType::permutation)
        {
            int numberOfShapes = 1;
            for (auto& entry : metadata)
            {
                numberOfShapes *= selectedMetadata(entry).size();
            }
            return numberOfShapes;
        }
        else if (rotamerType == GlycamMetadata::RotamerType::conformer)
        { // Conformer should mean that each dihedral will have the same number of metadata entries.
            // numberOfShapes = RotatableDihedrals_.size(); // This was correct for ASN for the wrong reason. 4
            // conformers and 4 dihedrals...
            return selectedMetadata(metadata[0]).size();
        }
        else
        {
            std::string str =
                "Error: Unknown rotamer type: " +
                std::to_string(static_cast<std::underlying_type_t<GlycamMetadata::RotamerType>>(rotamerType));
            gmml::log(__LINE__, __FILE__, gmml::ERR, str);
            throw std::runtime_error(str);
        }
    }
} // namespace

std::vector<size_t> cds::rotatableDihedralsWithMultipleRotamers(const std::vector<DihedralAngleDataVector>& metadata)
{
    std::vector<size_t> indices;
    for (size_t n = 0; n < metadata.size(); n++)
    {
        if (metadata[n].size() > 1)
        {
            indices.push_back(n);
        }
    }
    return indices;
}

size_t cds::numberOfShapes(GlycamMetadata::RotamerType rotamerType,
                           const std::vector<DihedralAngleDataVector>& metadata)
{
    auto selectMetadata = [](const DihedralAngleDataVector& data)
    {
        return data;
    };
    return getNumberOfShapesBase(selectMetadata, rotamerType, metadata);
}

size_t cds::numberOfLikelyShapes(GlycamMetadata::RotamerType rotamerType,
                                 const std::vector<DihedralAngleDataVector>& metadata)
{
    return getNumberOfShapesBase(GlycamMetadata::likelyMetadata, rotamerType, metadata);
}

cds::DihedralCoordinates cds::dihedralCoordinates(const cds::RotatableDihedral& dihedral)
{
    auto& atoms                       = dihedral.atoms;
    std::array<Coordinate*, 4> coords = {atoms[3]->getCoordinate(), atoms[2]->getCoordinate(),
                                         atoms[1]->getCoordinate(), atoms[0]->getCoordinate()};
    return DihedralCoordinates {*coords[0], *coords[1], *coords[2], *coords[3]};
}

std::string cds::print(const ResidueLink& link)
{
    std::stringstream ss;
    auto& linkageResidues = link.residues;
    ss << "ids: " << linkageResidues.first->getStringId() << "@" << link.atoms.first->getName() << " -- "
       << linkageResidues.second->getStringId() << "@" << link.atoms.second->getName() << "\n";
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

std::string cds::print(const ResidueLinkage& linkage)
{
    std::stringstream ss;
    ss << "ResidueLinkage Index: " << linkage.index << ", Name: " << linkage.name
       << ", NumberOfShapes: " << numberOfShapes(linkage.rotamerType, linkage.dihedralMetadata) << ", "
       << print(linkage.link);
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    for (auto& rotatableDihedral : linkage.rotatableDihedrals)
    {
        ss << print(rotatableDihedral);
    }
    return ss.str();
}
