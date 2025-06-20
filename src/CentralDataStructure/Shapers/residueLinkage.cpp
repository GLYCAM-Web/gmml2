#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <functional>
#include <sstream>

using cds::ResidueLinkage;
using cds::RotatableDihedral;

namespace
{
    int getNumberOfShapesBase(
        std::function<std::vector<size_t>(const GlycamMetadata::DihedralAngleDataTable&, const std::vector<size_t>&)>
            selectedMetadata,
        const GlycamMetadata::DihedralAngleDataTable& metadataTable, GlycamMetadata::RotamerType rotamerType,
        const std::vector<std::vector<size_t>>& metadata)
    {
        if (rotamerType == GlycamMetadata::RotamerType::permutation)
        {
            int numberOfShapes = 1;
            for (auto& entry : metadata)
            {
                numberOfShapes *= selectedMetadata(metadataTable, entry).size();
            }
            return numberOfShapes;
        }
        else if (rotamerType == GlycamMetadata::RotamerType::conformer)
        { // Conformer should mean that each dihedral will have the same number of metadata entries.
            // numberOfShapes = RotatableDihedrals_.size(); // This was correct for ASN for the wrong reason. 4
            // conformers and 4 dihedrals...
            return selectedMetadata(metadataTable, metadata[0]).size();
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

std::vector<size_t> cds::rotatableDihedralsWithMultipleRotamers(const std::vector<std::vector<size_t>>& metadata)
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

std::vector<cds::ResidueLinkage> cds::nonDerivativeResidueLinkages(const std::vector<cds::ResidueLinkage>& linkages)
{
    auto nonDerivative = [](const cds::ResidueLinkage& linkage)
    {
        return !linkage.isDerivative;
    };
    std::vector<cds::ResidueLinkage> result;
    result.reserve(linkages.size());
    std::copy_if(linkages.begin(), linkages.end(), std::back_inserter(result), nonDerivative);
    return result;
}

size_t cds::numberOfShapes(const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                           GlycamMetadata::RotamerType rotamerType, const std::vector<std::vector<size_t>>& metadata)
{
    auto selectMetadata = [](const GlycamMetadata::DihedralAngleDataTable&, const std::vector<size_t>& data)
    {
        return data;
    };
    return getNumberOfShapesBase(selectMetadata, metadataTable, rotamerType, metadata);
}

size_t cds::numberOfLikelyShapes(const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                 GlycamMetadata::RotamerType rotamerType,
                                 const std::vector<std::vector<size_t>>& metadata)
{
    return getNumberOfShapesBase(GlycamMetadata::likelyMetadata, metadataTable, rotamerType, metadata);
}

cds::DihedralCoordinates cds::dihedralCoordinates(const cds::RotatableDihedral& dihedral)
{
    auto& atoms = dihedral.atoms;
    return {atoms[3]->coordinate(), atoms[2]->coordinate(), atoms[1]->coordinate(), atoms[0]->coordinate()};
}

std::string cds::print(const ResidueLink& link)
{
    std::stringstream ss;
    auto& linkageResidues = link.residues;
    ss << "ids: " << cds::residueStringId(linkageResidues.first) << "@" << link.atoms.first->getName() << " -- "
       << cds::residueStringId(linkageResidues.second) << "@" << link.atoms.second->getName() << "\n";
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

std::string cds::print(const GlycamMetadata::DihedralAngleDataTable& table, const ResidueLinkage& linkage)
{
    std::stringstream ss;
    ss << "ResidueLinkage Index: " << linkage.index << ", Name: " << linkage.name
       << ", NumberOfShapes: " << numberOfShapes(table, linkage.rotamerType, linkage.dihedralMetadata) << ", "
       << print(linkage.link);
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    for (auto& rotatableDihedral : linkage.rotatableDihedrals)
    {
        ss << print(rotatableDihedral);
    }
    return ss.str();
}
