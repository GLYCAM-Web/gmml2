#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomCoordinates.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <variant>

using cds::RotatableDihedral;

void cds::setDihedralAngle(RotatableDihedral& dihedral, cds::AngleWithMetadata target)
{
    auto matrix = rotationTo(dihedralCoordinates(dihedral), constants::toRadians(target.value));
    matrix.rotateCoordinates(atomCoordinateReferences(dihedral.movingAtoms));
    dihedral.currentMetadataIndex = target.metadataIndex;
}

bool cds::setSpecificShape(RotatableDihedral& dihedral, const DihedralAngleDataVector& metadataVector,
                           std::string dihedralName, std::string selectedRotamer)
{
    if (dihedralName == metadataVector[0].dihedral_angle_name_)
    {
        for (size_t n = 0; n < metadataVector.size(); n++)
        {
            auto& metadata = metadataVector[n];
            if (metadata.rotamer_name_ == selectedRotamer)
            {
                setDihedralAngle(dihedral, {metadata.default_angle_value_, metadata.default_angle_value_, n});
                return true;
            }
        }
    }
    return false;
}

void cds::setSpecificShape(std::vector<RotatableDihedral>& dihedrals,
                           const std::vector<DihedralAngleDataVector>& metadata, std::string dihedralName,
                           std::string selectedRotamer)
{
    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        // This will call RotatableDihedrals that don't have dihedralName (phi,psi), and nothing will happen. Hmmm.
        if (setSpecificShape(dihedrals[n], metadata[n], dihedralName, selectedRotamer))
        {
            return; // Return once you manage to set a shape.
        }
    }
    std::string errorMessage = "Did not set " + dihedralName + " to " + selectedRotamer +
                               " as requested in ResidueLinkage::SetSpecificShape()";
    gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
    throw std::runtime_error(errorMessage);
}

std::vector<cds::AngleWithMetadata> cds::currentShape(const std::vector<RotatableDihedral>& dihedrals,
                                                      const std::vector<DihedralAngleDataVector>& metadata)
{
    std::vector<AngleWithMetadata> result;
    result.reserve(dihedrals.size());

    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        auto& dihedral        = dihedrals[n];
        auto& metadataIndex   = dihedral.currentMetadataIndex;
        auto& currentMetadata = metadata[n][metadataIndex];
        result.push_back({constants::toDegrees(cds::angle(dihedralCoordinates(dihedral))),
                          currentMetadata.default_angle_value_, metadataIndex});
    }

    return result;
}

std::vector<std::vector<cds::AngleWithMetadata>> cds::currentShape(const std::vector<ResidueLinkage>& linkages)
{
    std::vector<std::vector<AngleWithMetadata>> result;
    result.reserve(linkages.size());

    for (auto& a : linkages)
    {
        result.push_back(currentShape(a.rotatableDihedrals, a.dihedralMetadata));
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
    auto& dihedrals = linkage.rotatableDihedrals;
    if (linkage.rotamerType == GlycamMetadata::RotamerType::conformer)
    {
        if (!std::holds_alternative<ConformerShapePreference>(preference))
        {
            throw std::runtime_error("expected but did not receive ConformerShapePreference for conformer linkage");
        }
        ConformerShapePreference pref = std::get<ConformerShapePreference>(preference);
        size_t metadataIndex          = pref.metadataOrder[0];
        for (size_t n = 0; n < dihedrals.size(); n++)
        {
            double angle = pref.angles[n][metadataIndex];
            setDihedralAngle(dihedrals[n], {angle, angle, metadataIndex});
        }
    }
    else
    {
        if (!std::holds_alternative<PermutationShapePreference>(preference))
        {
            throw std::runtime_error("expected but did not receive PermutationShapePreference for permutation linkage");
        }
        PermutationShapePreference pref = std::get<PermutationShapePreference>(preference);
        for (size_t n = 0; n < dihedrals.size(); n++)
        {
            size_t metadataIndex = pref.metadataOrder[n][0];
            double angle         = pref.angles[n][metadataIndex];
            setDihedralAngle(dihedrals[n], {angle, angle, metadataIndex});
        }
    }
}

void cds::setShapeToPreference(std::vector<ResidueLinkage>& linkages,
                               const std::vector<ResidueLinkageShapePreference>& preferences)
{
    for (size_t n = 0; n < linkages.size(); n++)
    {
        setShapeToPreference(linkages[n], preferences[n]);
    }
}

cds::ResidueLinkageShapePreference cds::linkageShapePreference(MetadataDistribution metadataDistribution,
                                                               AngleDistribution angleDistribution,
                                                               const ResidueLinkage& linkage)
{
    std::vector<std::vector<double>> angles;
    angles.resize(linkage.dihedralMetadata.size());
    for (size_t n = 0; n < linkage.dihedralMetadata.size(); n++)
    {
        for (auto& metadata : linkage.dihedralMetadata[n])
        {
            angles[n].push_back(angleDistribution(metadata));
        }
    }
    if (linkage.rotamerType == GlycamMetadata::RotamerType::conformer)
    {
        auto order    = metadataDistribution(linkage.dihedralMetadata[0]);
        auto isFrozen = std::vector<bool>(linkage.rotatableDihedrals.size(), false);
        return ConformerShapePreference {isFrozen, angles, order};
    }
    else
    {
        std::vector<std::vector<size_t>> order;
        order.reserve(linkage.dihedralMetadata.size());
        for (auto& metadataVector : linkage.dihedralMetadata)
        {
            order.push_back(metadataDistribution(metadataVector));
        }
        return PermutationShapePreference {angles, order};
    }
}

cds::ResidueLinkageShapePreference cds::defaultShapePreference(const ResidueLinkage& linkage)
{
    auto metadataOrder = [](const DihedralAngleDataVector metadataVector)
    {
        return codeUtils::indexVector(metadataVector);
    };
    auto defaultAngle = [](const DihedralAngleData metadata)
    {
        return metadata.default_angle_value_;
    };
    return linkageShapePreference(metadataOrder, defaultAngle, linkage);
}

cds::ResidueLinkageShapePreference cds::selectedRotamersOnly(MetadataPreferenceSelection metadataSelection,
                                                             const ResidueLinkage& linkage,
                                                             const ResidueLinkageShapePreference& preference)
{
    if (std::holds_alternative<ConformerShapePreference>(preference))
    {
        auto pref     = std::get<ConformerShapePreference>(preference);
        auto isFrozen = std::vector<bool>(linkage.rotatableDihedrals.size(), false);
        return ConformerShapePreference {isFrozen, pref.angles,
                                         metadataSelection(pref.metadataOrder, linkage.rotatableDihedrals[0])};
    }
    else if (std::holds_alternative<PermutationShapePreference>(preference))
    {
        auto pref  = std::get<PermutationShapePreference>(preference);
        auto order = pref.metadataOrder;
        std::vector<std::vector<size_t>> selected;
        selected.reserve(order.size());
        for (size_t n = 0; n < order.size(); n++)
        {
            selected.push_back(metadataSelection(order[n], linkage.rotatableDihedrals[n]));
        }
        return PermutationShapePreference {pref.angles, selected};
    }
    else
    {
        throw std::runtime_error("unhandled linkage shape preference in cds::selectedRotamersOnly");
    }
}

cds::ResidueLinkageShapePreference cds::firstRotamerOnly(const ResidueLinkage& linkage,
                                                         const ResidueLinkageShapePreference& preference)
{
    auto first = [](std::vector<size_t> order, const RotatableDihedral&)
    {
        return std::vector<size_t> {order[0]};
    };
    return selectedRotamersOnly(first, linkage, preference);
}

cds::ResidueLinkageShapePreference cds::currentRotamerOnly(const ResidueLinkage& linkage,
                                                           const ResidueLinkageShapePreference& preference)
{
    auto current = [](std::vector<size_t>, const RotatableDihedral& dihedral)
    {
        return std::vector<size_t> {dihedral.currentMetadataIndex};
    };
    return selectedRotamersOnly(current, linkage, preference);
}

std::vector<cds::ResidueLinkageShapePreference> cds::linkageShapePreference(MetadataDistribution metadataDistribution,
                                                                            AngleDistribution angleDistribution,
                                                                            const std::vector<ResidueLinkage>& linkages)
{
    std::vector<ResidueLinkageShapePreference> result;
    result.reserve(linkages.size());
    for (auto& linkage : linkages)
    {
        result.push_back(linkageShapePreference(metadataDistribution, angleDistribution, linkage));
    }
    return result;
}

std::vector<cds::ResidueLinkageShapePreference> cds::defaultShapePreference(const std::vector<ResidueLinkage>& linkages)
{
    std::vector<ResidueLinkageShapePreference> result;
    result.reserve(linkages.size());
    for (auto& linkage : linkages)
    {
        result.push_back(defaultShapePreference(linkage));
    }
    return result;
}

std::vector<cds::ResidueLinkageShapePreference>
cds::firstRotamerOnly(const std::vector<ResidueLinkage>& linkages,
                      const std::vector<ResidueLinkageShapePreference>& preferences)
{
    std::vector<ResidueLinkageShapePreference> result;
    result.reserve(preferences.size());
    for (size_t n = 0; n < preferences.size(); n++)
    {
        result.push_back(firstRotamerOnly(linkages[n], preferences[n]));
    }
    return result;
}

std::vector<cds::ResidueLinkageShapePreference>
cds::currentRotamerOnly(const std::vector<ResidueLinkage>& linkages,
                        const std::vector<ResidueLinkageShapePreference>& preferences)
{
    std::vector<ResidueLinkageShapePreference> result;
    result.reserve(preferences.size());
    for (size_t n = 0; n < preferences.size(); n++)
    {
        result.push_back(currentRotamerOnly(linkages[n], preferences[n]));
    }
    return result;
}
