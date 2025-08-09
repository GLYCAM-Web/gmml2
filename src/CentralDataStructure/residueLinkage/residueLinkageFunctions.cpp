#include "include/CentralDataStructure/residueLinkage/residueLinkageFunctions.hpp"

#include "include/CentralDataStructure/atom.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/orientation.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/util/logging.hpp"

#include <functional>
#include <sstream>

namespace gmml
{
    namespace
    {
        int getNumberOfShapesBase(
            std::function<std::vector<size_t>(const DihedralAngleDataTable&, const std::vector<size_t>&)>
                selectedMetadata,
            const DihedralAngleDataTable& metadataTable,
            RotamerType rotamerType,
            const std::vector<std::vector<size_t>>& metadata)
        {
            if (rotamerType == RotamerType::permutation)
            {
                int numberOfShapes = 1;
                for (auto& entry : metadata)
                {
                    numberOfShapes *= selectedMetadata(metadataTable, entry).size();
                }
                return numberOfShapes;
            }
            else if (rotamerType == RotamerType::conformer)
            { // Conformer should mean that each dihedral will have the same number of metadata entries.
                // numberOfShapes = RotatableDihedrals_.size(); // This was correct for ASN for the wrong reason. 4
                // conformers and 4 dihedrals...
                return selectedMetadata(metadataTable, metadata[0]).size();
            }
            else
            {
                std::string str = "Error: Unknown rotamer type: " +
                                  std::to_string(static_cast<std::underlying_type_t<RotamerType>>(rotamerType));
                util::log(__LINE__, __FILE__, util::ERR, str);
                throw std::runtime_error(str);
            }
        }
    } // namespace

    std::vector<size_t> rotatableDihedralsWithMultipleRotamers(const std::vector<std::vector<size_t>>& metadata)
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

    std::vector<ResidueLinkage> nonDerivativeResidueLinkages(const std::vector<ResidueLinkage>& linkages)
    {
        auto nonDerivative = [](const ResidueLinkage& linkage) { return !linkage.isDerivative; };
        std::vector<ResidueLinkage> result;
        result.reserve(linkages.size());
        std::copy_if(linkages.begin(), linkages.end(), std::back_inserter(result), nonDerivative);
        return result;
    }

    size_t numberOfShapes(
        const DihedralAngleDataTable& metadataTable,
        RotamerType rotamerType,
        const std::vector<std::vector<size_t>>& metadata)
    {
        auto selectMetadata = [](const DihedralAngleDataTable&, const std::vector<size_t>& data) { return data; };
        return getNumberOfShapesBase(selectMetadata, metadataTable, rotamerType, metadata);
    }

    size_t numberOfLikelyShapes(
        const DihedralAngleDataTable& metadataTable,
        RotamerType rotamerType,
        const std::vector<std::vector<size_t>>& metadata)
    {
        return getNumberOfShapesBase(likelyMetadata, metadataTable, rotamerType, metadata);
    }

    DihedralCoordinates dihedralCoordinates(const RotatableBond& bond)
    {
        auto& atoms = bond.dihedralAtoms;
        return {atoms[3]->coordinate(), atoms[2]->coordinate(), atoms[1]->coordinate(), atoms[0]->coordinate()};
    }

    std::string print(const ResidueLink& link)
    {
        std::stringstream ss;
        auto& linkageResidues = link.residues;
        ss << "ids: " << residueStringId(linkageResidues.first) << "@" << link.atoms.first->getName() << " -- "
           << residueStringId(linkageResidues.second) << "@" << link.atoms.second->getName() << "\n";
        return ss.str();
    }

    std::string print(const RotatableBond& bond)
    {
        auto& atoms = bond.dihedralAtoms;
        std::stringstream ss;
        ss << atoms[0]->getName() << ", " << atoms[1]->getName() << ", " << atoms[2]->getName() << ", "
           << atoms[3]->getName() << ": " << angle(dihedralCoordinates(bond)) << ".\n";
        util::log(__LINE__, __FILE__, util::INF, ss.str());
        return ss.str();
    }

    std::string print(const DihedralAngleDataTable& table, const ResidueLinkage& linkage)
    {
        std::stringstream ss;
        ss << "ResidueLinkage Index: " << linkage.index << ", Name: " << linkage.name
           << ", NumberOfShapes: " << numberOfShapes(table, linkage.rotamerType, linkage.dihedralMetadata) << ", "
           << print(linkage.link);
        util::log(__LINE__, __FILE__, util::INF, ss.str());
        for (auto& rotatableDihedral : linkage.rotatableBonds)
        {
            ss << print(rotatableDihedral);
        }
        return ss.str();
    }
} // namespace gmml
