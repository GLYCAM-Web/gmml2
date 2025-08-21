#include "include/carbohydrate/dihedralShape.hpp"

#include "include/CentralDataStructure/atom.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/residueLinkage/residueLinkageFunctions.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/matrix.hpp"
#include "include/util/constants.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"

#include <string>
#include <variant>
#include <vector>

namespace gmml
{
    void setDihedralAngle(RotatableBond& bond, AngleWithMetadata target)
    {
        Matrix4x4 matrix = rotationTo(dihedralCoordinates(bond), constants::toRadians(target.value));
        setAtomCoordinates(bond.movingAtoms, transform(matrix, atomCoordinates(bond.movingAtoms)));
        bond.currentMetadataIndex = target.metadataIndex;
    }

    bool setSpecificShape(
        const DihedralAngleDataTable& metadataTable,
        RotatableBond& bond,
        const std::vector<size_t>& metadataVector,
        std::string dihedralName,
        std::string selectedRotamer)
    {
        if (dihedralName == metadataTable.entries[metadataVector[0]].dihedral_angle_name_)
        {
            for (size_t n = 0; n < metadataVector.size(); n++)
            {
                auto& metadata = metadataTable.entries[metadataVector[n]];
                if (metadata.rotamer_name_ == selectedRotamer)
                {
                    setDihedralAngle(bond, {metadata.default_angle, metadata.default_angle, n});
                    return true;
                }
            }
        }
        return false;
    }

    void setSpecificShape(
        const DihedralAngleDataTable& metadataTable,
        std::vector<RotatableBond>& bonds,
        const std::vector<std::vector<size_t>>& metadata,
        std::string dihedralName,
        std::string selectedRotamer)
    {
        for (size_t n = 0; n < bonds.size(); n++)
        {
            // This will call RotatableDihedrals that don't have dihedralName (phi,psi), and nothing will happen. Hmmm.
            if (setSpecificShape(metadataTable, bonds[n], metadata[n], dihedralName, selectedRotamer))
            {
                return; // Return once you manage to set a shape.
            }
        }
        std::string errorMessage = "Did not set " + dihedralName + " to " + selectedRotamer +
                                   " as requested in ResidueLinkage::SetSpecificShape()";
        util::log(__LINE__, __FILE__, util::ERR, errorMessage);
        throw std::runtime_error(errorMessage);
    }

    std::vector<AngleWithMetadata> currentShape(
        const DihedralAngleDataTable& metadataTable,
        const std::vector<RotatableBond>& bonds,
        const std::vector<std::vector<size_t>>& metadata)
    {
        std::vector<AngleWithMetadata> result;
        result.reserve(bonds.size());

        for (size_t n = 0; n < bonds.size(); n++)
        {
            auto& bond = bonds[n];
            auto& metadataIndex = bond.currentMetadataIndex;
            auto& currentMetadata = metadata[n][metadataIndex];
            result.push_back(
                {constants::toDegrees(angle(dihedralCoordinates(bond))),
                 metadataTable.entries[currentMetadata].default_angle,
                 metadataIndex});
        }

        return result;
    }

    std::vector<std::vector<AngleWithMetadata>> currentShape(
        const DihedralAngleDataTable& metadataTable, const std::vector<ResidueLinkage>& linkages)
    {
        std::vector<std::vector<AngleWithMetadata>> result;
        result.reserve(linkages.size());

        for (auto& a : linkages)
        {
            result.push_back(currentShape(metadataTable, a.rotatableBonds, a.dihedralMetadata));
        }

        return result;
    }

    void setShape(std::vector<ResidueLinkage>& linkages, const std::vector<std::vector<AngleWithMetadata>>& angles)
    {
        for (size_t n = 0; n < linkages.size(); n++)
        {
            setShape(linkages[n].rotatableBonds, angles[n]);
        }
    }

    void setShape(std::vector<RotatableBond>& bonds, const std::vector<AngleWithMetadata>& angles)
    {
        for (size_t n = 0; n < angles.size(); n++)
        {
            setDihedralAngle(bonds[n], angles[n]);
        }
    }

    void setShapeToPreference(ResidueLinkage& linkage, const ResidueLinkageShapePreference& preference)
    {
        std::vector<RotatableBond>& bonds = linkage.rotatableBonds;
        std::function<void(const ConformerShapePreference&)> onConformer = [&](const ConformerShapePreference& pref)
        {
            if (linkage.rotamerType != RotamerType::conformer)
            {
                throw std::runtime_error("expected but did not receive ConformerShapePreference for conformer linkage");
            }
            size_t metadataIndex = pref.metadataOrder[0];
            for (size_t n = 0; n < bonds.size(); n++)
            {
                double angle = pref.angles[n][metadataIndex];
                setDihedralAngle(bonds[n], {angle, angle, metadataIndex});
            }
        };
        std::function<void(const PermutationShapePreference&)> onPermutation =
            [&](const PermutationShapePreference& pref)
        {
            if (linkage.rotamerType != RotamerType::permutation)
            {
                throw std::runtime_error(
                    "expected but did not receive PermutationShapePreference for permutation linkage");
            }
            for (size_t n = 0; n < bonds.size(); n++)
            {
                size_t metadataIndex = pref.metadataOrder[n][0];
                double angle = pref.angles[n][metadataIndex];
                setDihedralAngle(bonds[n], {angle, angle, metadataIndex});
            }
        };
        return onResidueLinkageShapePreference(onConformer, onPermutation, preference);
    }

    void setShapeToPreference(std::vector<ResidueLinkage>& linkages, const GlycanShapePreference& preferences)
    {
        for (size_t n = 0; n < linkages.size(); n++)
        {
            setShapeToPreference(linkages[n], preferences[n]);
        }
    }

    ResidueLinkageShapePreference linkageShapePreference(
        MetadataDistribution metadataDistribution,
        AngleDistribution angleDistribution,
        const DihedralAngleDataTable& metadataTable,
        const RotamerType rotamerType,
        const std::vector<std::vector<size_t>>& dihedralMetadata)
    {
        std::vector<std::vector<double>> angles;
        angles.resize(dihedralMetadata.size());
        for (size_t n = 0; n < dihedralMetadata.size(); n++)
        {
            for (auto& metadata : dihedralMetadata[n])
            {
                angles[n].push_back(angleDistribution(metadataTable.entries[metadata]));
            }
        }
        if (rotamerType == RotamerType::conformer)
        {
            auto order = metadataDistribution(metadataTable, dihedralMetadata[0]);
            auto isFrozen = std::vector<bool>(dihedralMetadata.size(), false);
            return ConformerShapePreference {isFrozen, angles, order};
        }
        else
        {
            std::vector<std::vector<size_t>> order;
            order.reserve(dihedralMetadata.size());
            for (auto& metadataVector : dihedralMetadata)
            {
                order.push_back(metadataDistribution(metadataTable, metadataVector));
            }
            return PermutationShapePreference {angles, order};
        }
    }

    ResidueLinkageShapePreference defaultShapePreference(
        const DihedralAngleDataTable& metadataTable,
        const RotamerType rotamerType,
        const std::vector<std::vector<size_t>>& dihedralMetadata)
    {
        auto metadataOrder = [](const DihedralAngleDataTable&, const std::vector<size_t> metadataVector)
        { return util::indexVector(metadataVector); };
        auto defaultAngle = [](const DihedralAngleData& metadata) { return metadata.default_angle; };
        return linkageShapePreference(metadataOrder, defaultAngle, metadataTable, rotamerType, dihedralMetadata);
    }

    ResidueLinkageShapePreference selectedRotamersOnly(
        MetadataPreferenceSelection metadataSelection,
        const ResidueLinkage& linkage,
        const ResidueLinkageShapePreference& preference)
    {
        std::function<ResidueLinkageShapePreference(const ConformerShapePreference&)> onConformer =
            [&](const ConformerShapePreference& pref)
        {
            auto isFrozen = std::vector<bool>(linkage.rotatableBonds.size(), false);
            return ConformerShapePreference {
                isFrozen, pref.angles, metadataSelection(pref.metadataOrder, linkage.rotatableBonds[0])};
        };
        std::function<ResidueLinkageShapePreference(const PermutationShapePreference&)> onPermutation =
            [&](const PermutationShapePreference& pref)
        {
            auto order = pref.metadataOrder;
            std::vector<std::vector<size_t>> selected;
            selected.reserve(order.size());
            for (size_t n = 0; n < order.size(); n++)
            {
                selected.push_back(metadataSelection(order[n], linkage.rotatableBonds[n]));
            }
            return PermutationShapePreference {pref.angles, selected};
        };
        return onResidueLinkageShapePreference(onConformer, onPermutation, preference);
    }

    ResidueLinkageShapePreference firstRotamerOnly(
        const ResidueLinkage& linkage, const ResidueLinkageShapePreference& preference)
    {
        auto first = [](std::vector<size_t> order, const RotatableBond&) { return std::vector<size_t> {order[0]}; };
        return selectedRotamersOnly(first, linkage, preference);
    }

    ResidueLinkageShapePreference currentRotamerOnly(
        const ResidueLinkage& linkage, const ResidueLinkageShapePreference& preference)
    {
        auto current = [](std::vector<size_t>, const RotatableBond& bond)
        { return std::vector<size_t> {bond.currentMetadataIndex}; };
        return selectedRotamersOnly(current, linkage, preference);
    }

    GlycanShapePreference linkageShapePreference(
        MetadataDistribution metadataDistribution,
        AngleDistribution angleDistribution,
        const DihedralAngleDataTable& metadataTable,
        const std::vector<ResidueLinkage>& linkages)
    {
        GlycanShapePreference result;
        result.reserve(linkages.size());
        for (auto& linkage : linkages)
        {
            result.push_back(linkageShapePreference(
                metadataDistribution, angleDistribution, metadataTable, linkage.rotamerType, linkage.dihedralMetadata));
        }
        return result;
    }

    GlycanShapePreference defaultShapePreference(
        const DihedralAngleDataTable& metadataTable, const std::vector<ResidueLinkage>& linkages)
    {
        GlycanShapePreference result;
        result.reserve(linkages.size());
        for (auto& linkage : linkages)
        {
            result.push_back(defaultShapePreference(metadataTable, linkage.rotamerType, linkage.dihedralMetadata));
        }
        return result;
    }

    GlycanShapePreference firstRotamerOnly(
        const std::vector<ResidueLinkage>& linkages, const GlycanShapePreference& preferences)
    {
        GlycanShapePreference result;
        result.reserve(preferences.size());
        for (size_t n = 0; n < preferences.size(); n++)
        {
            result.push_back(firstRotamerOnly(linkages[n], preferences[n]));
        }
        return result;
    }

    GlycanShapePreference currentRotamerOnly(
        const std::vector<ResidueLinkage>& linkages, const GlycanShapePreference& preferences)
    {
        GlycanShapePreference result;
        result.reserve(preferences.size());
        for (size_t n = 0; n < preferences.size(); n++)
        {
            result.push_back(currentRotamerOnly(linkages[n], preferences[n]));
        }
        return result;
    }
} // namespace gmml
