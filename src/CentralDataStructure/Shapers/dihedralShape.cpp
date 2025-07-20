#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageFunctions.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <variant>

using cds::RotatableDihedral;

void cds::setDihedralAngle(RotatableDihedral& dihedral, cds::AngleWithMetadata target)
{
    cds::RotationMatrix matrix = rotationTo(dihedralCoordinates(dihedral), constants::toRadians(target.value));
    std::vector<Coordinate> coordinates = atomCoordinates(dihedral.movingAtoms);
    matrix.rotateCoordinates(coordinates);
    setAtomCoordinates(dihedral.movingAtoms, coordinates);
    dihedral.currentMetadataIndex = target.metadataIndex;
}

bool cds::setSpecificShape(
    const GlycamMetadata::DihedralAngleDataTable& metadataTable,
    RotatableDihedral& dihedral,
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
                setDihedralAngle(dihedral, {metadata.default_angle, metadata.default_angle, n});
                return true;
            }
        }
    }
    return false;
}

void cds::setSpecificShape(
    const GlycamMetadata::DihedralAngleDataTable& metadataTable,
    std::vector<RotatableDihedral>& dihedrals,
    const std::vector<std::vector<size_t>>& metadata,
    std::string dihedralName,
    std::string selectedRotamer)
{
    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        // This will call RotatableDihedrals that don't have dihedralName (phi,psi), and nothing will happen. Hmmm.
        if (setSpecificShape(metadataTable, dihedrals[n], metadata[n], dihedralName, selectedRotamer))
        {
            return; // Return once you manage to set a shape.
        }
    }
    std::string errorMessage = "Did not set " + dihedralName + " to " + selectedRotamer +
                               " as requested in ResidueLinkage::SetSpecificShape()";
    gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
    throw std::runtime_error(errorMessage);
}

std::vector<cds::AngleWithMetadata> cds::currentShape(
    const GlycamMetadata::DihedralAngleDataTable& metadataTable,
    const std::vector<RotatableDihedral>& dihedrals,
    const std::vector<std::vector<size_t>>& metadata)
{
    std::vector<AngleWithMetadata> result;
    result.reserve(dihedrals.size());

    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        auto& dihedral = dihedrals[n];
        auto& metadataIndex = dihedral.currentMetadataIndex;
        auto& currentMetadata = metadata[n][metadataIndex];
        result.push_back(
            {constants::toDegrees(cds::angle(dihedralCoordinates(dihedral))),
             metadataTable.entries[currentMetadata].default_angle,
             metadataIndex});
    }

    return result;
}

std::vector<std::vector<cds::AngleWithMetadata>> cds::currentShape(
    const GlycamMetadata::DihedralAngleDataTable& metadataTable, const std::vector<ResidueLinkage>& linkages)
{
    std::vector<std::vector<AngleWithMetadata>> result;
    result.reserve(linkages.size());

    for (auto& a : linkages)
    {
        result.push_back(currentShape(metadataTable, a.rotatableDihedrals, a.dihedralMetadata));
    }

    return result;
}

void cds::setShape(std::vector<ResidueLinkage>& linkages, const std::vector<std::vector<AngleWithMetadata>>& angles)
{
    for (size_t n = 0; n < linkages.size(); n++)
    {
        setShape(linkages[n].rotatableDihedrals, angles[n]);
    }
}

void cds::setShape(std::vector<cds::RotatableDihedral>& dihedrals, const std::vector<AngleWithMetadata>& angles)
{
    for (size_t n = 0; n < angles.size(); n++)
    {
        setDihedralAngle(dihedrals[n], angles[n]);
    }
}

void cds::setShapeToPreference(ResidueLinkage& linkage, const ResidueLinkageShapePreference& preference)
{
    std::vector<RotatableDihedral>& dihedrals = linkage.rotatableDihedrals;
    std::function<void(const ConformerShapePreference&)> onConformer = [&](const ConformerShapePreference& pref)
    {
        if (linkage.rotamerType != GlycamMetadata::RotamerType::conformer)
        {
            throw std::runtime_error("expected but did not receive ConformerShapePreference for conformer linkage");
        }
        size_t metadataIndex = pref.metadataOrder[0];
        for (size_t n = 0; n < dihedrals.size(); n++)
        {
            double angle = pref.angles[n][metadataIndex];
            setDihedralAngle(dihedrals[n], {angle, angle, metadataIndex});
        }
    };
    std::function<void(const PermutationShapePreference&)> onPermutation = [&](const PermutationShapePreference& pref)
    {
        if (linkage.rotamerType != GlycamMetadata::RotamerType::permutation)
        {
            throw std::runtime_error("expected but did not receive PermutationShapePreference for permutation linkage");
        }
        for (size_t n = 0; n < dihedrals.size(); n++)
        {
            size_t metadataIndex = pref.metadataOrder[n][0];
            double angle = pref.angles[n][metadataIndex];
            setDihedralAngle(dihedrals[n], {angle, angle, metadataIndex});
        }
    };
    return onResidueLinkageShapePreference(onConformer, onPermutation, preference);
}

void cds::setShapeToPreference(std::vector<ResidueLinkage>& linkages, const GlycanShapePreference& preferences)
{
    for (size_t n = 0; n < linkages.size(); n++)
    {
        setShapeToPreference(linkages[n], preferences[n]);
    }
}

cds::ResidueLinkageShapePreference cds::linkageShapePreference(
    MetadataDistribution metadataDistribution,
    AngleDistribution angleDistribution,
    const GlycamMetadata::DihedralAngleDataTable& metadataTable,
    const GlycamMetadata::RotamerType rotamerType,
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
    if (rotamerType == GlycamMetadata::RotamerType::conformer)
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

cds::ResidueLinkageShapePreference cds::defaultShapePreference(
    const GlycamMetadata::DihedralAngleDataTable& metadataTable,
    const GlycamMetadata::RotamerType rotamerType,
    const std::vector<std::vector<size_t>>& dihedralMetadata)
{
    auto metadataOrder = [](const GlycamMetadata::DihedralAngleDataTable&, const std::vector<size_t> metadataVector)
    { return codeUtils::indexVector(metadataVector); };
    auto defaultAngle = [](const GlycamMetadata::DihedralAngleData& metadata) { return metadata.default_angle; };
    return linkageShapePreference(metadataOrder, defaultAngle, metadataTable, rotamerType, dihedralMetadata);
}

cds::ResidueLinkageShapePreference cds::selectedRotamersOnly(
    MetadataPreferenceSelection metadataSelection,
    const ResidueLinkage& linkage,
    const ResidueLinkageShapePreference& preference)
{
    std::function<ResidueLinkageShapePreference(const ConformerShapePreference&)> onConformer =
        [&](const ConformerShapePreference& pref)
    {
        auto isFrozen = std::vector<bool>(linkage.rotatableDihedrals.size(), false);
        return ConformerShapePreference {
            isFrozen, pref.angles, metadataSelection(pref.metadataOrder, linkage.rotatableDihedrals[0])};
    };
    std::function<ResidueLinkageShapePreference(const PermutationShapePreference&)> onPermutation =
        [&](const PermutationShapePreference& pref)
    {
        auto order = pref.metadataOrder;
        std::vector<std::vector<size_t>> selected;
        selected.reserve(order.size());
        for (size_t n = 0; n < order.size(); n++)
        {
            selected.push_back(metadataSelection(order[n], linkage.rotatableDihedrals[n]));
        }
        return PermutationShapePreference {pref.angles, selected};
    };
    return onResidueLinkageShapePreference(onConformer, onPermutation, preference);
}

cds::ResidueLinkageShapePreference cds::firstRotamerOnly(
    const ResidueLinkage& linkage, const ResidueLinkageShapePreference& preference)
{
    auto first = [](std::vector<size_t> order, const RotatableDihedral&) { return std::vector<size_t> {order[0]}; };
    return selectedRotamersOnly(first, linkage, preference);
}

cds::ResidueLinkageShapePreference cds::currentRotamerOnly(
    const ResidueLinkage& linkage, const ResidueLinkageShapePreference& preference)
{
    auto current = [](std::vector<size_t>, const RotatableDihedral& dihedral)
    { return std::vector<size_t> {dihedral.currentMetadataIndex}; };
    return selectedRotamersOnly(current, linkage, preference);
}

cds::GlycanShapePreference cds::linkageShapePreference(
    MetadataDistribution metadataDistribution,
    AngleDistribution angleDistribution,
    const GlycamMetadata::DihedralAngleDataTable& metadataTable,
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

cds::GlycanShapePreference cds::defaultShapePreference(
    const GlycamMetadata::DihedralAngleDataTable& metadataTable, const std::vector<ResidueLinkage>& linkages)
{
    GlycanShapePreference result;
    result.reserve(linkages.size());
    for (auto& linkage : linkages)
    {
        result.push_back(defaultShapePreference(metadataTable, linkage.rotamerType, linkage.dihedralMetadata));
    }
    return result;
}

cds::GlycanShapePreference cds::firstRotamerOnly(
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

cds::GlycanShapePreference cds::currentRotamerOnly(
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
