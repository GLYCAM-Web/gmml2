#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <functional>

namespace glycoproteinBuilder
{
    void updateGlycanBounds(const assembly::Graph& graph, const AssemblyData& data, assembly::Bounds& bounds,
                            size_t glycanId)
    {
        size_t siteResidue = data.glycans.attachmentResidue[glycanId];
        updateResidueBounds(graph, bounds, siteResidue);
        updateResidueMoleculeBounds(graph, bounds, siteResidue);
        size_t moleculeId = data.glycans.moleculeId[glycanId];
        for (size_t residue : moleculeResidues(graph, moleculeId))
        {
            updateResidueBounds(graph, bounds, residue);
        }
        updateMoleculeBounds(graph, bounds, moleculeId);
    }

    std::array<cds::Coordinate, 4> dihedralCoordinates(const AssemblyData& data, const assembly::Bounds& bounds,
                                                       size_t dihedralId)
    {
        const std::array<size_t, 4>& atoms = data.indices.rotatableDihedrals[dihedralId].atoms;
        auto coord                         = [&](size_t n)
        {
            return bounds.atoms[atoms[n]].center;
        };
        return {coord(3), coord(2), coord(1), coord(0)};
    }

    void setDihedralAngle(const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                          size_t dihedralId, const cds::AngleWithMetadata& target)
    {
        mutableData.dihedralCurrentMetadata[dihedralId] = target.metadataIndex;
        const RotatableDihedralIndices& dihedral        = data.indices.rotatableDihedrals[dihedralId];
        std::array<cds::Coordinate, 4> dihedralCoords   = dihedralCoordinates(data, mutableData.bounds, dihedralId);
        auto matrix = cds::rotationTo(dihedralCoords, constants::toRadians(target.value));
        for (size_t atom : dihedral.movingAtoms)
        {
            cds::Coordinate& coord = mutableData.bounds.atoms[atom].center;
            coord                  = matrix * coord;
        }
        updateBoundsContainingAtoms(graph, mutableData.bounds, dihedral.movingAtoms);
    }

    void setLinkageShape(const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                         size_t glycanId, const std::vector<cds::AngleWithMetadata>& recordedShape)
    {
        const std::vector<size_t>& linkages = data.glycans.linkages[glycanId];
        for (size_t n = 0; n < linkages.size(); n++)
        {
            size_t linkageId                     = linkages[n];
            const std::vector<size_t>& dihedrals = data.indices.residueLinkages[linkageId].rotatableDihedrals;
            for (size_t k = 0; k < dihedrals.size(); k++)
            {
                size_t dihedralId = dihedrals[k];
                setDihedralAngle(graph, data, mutableData, dihedralId, recordedShape[dihedralId]);
            }
        }
    }

    void setLinkageShapeToPreference(const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                                     size_t linkageId, const cds::ResidueLinkageShapePreference& preference)
    {
        GlycamMetadata::RotamerType rotamerType = data.residueLinkageData.rotamerTypes[linkageId];
        const std::vector<size_t>& dihedrals    = data.indices.residueLinkages[linkageId].rotatableDihedrals;

        std::function<void(const cds::ConformerShapePreference&)> onConformer =
            [&](const cds::ConformerShapePreference& pref)
        {
            if (rotamerType != GlycamMetadata::RotamerType::conformer)
            {
                throw std::runtime_error("expected but did not receive ConformerShapePreference for conformer linkage");
            }
            size_t metadataIndex = pref.metadataOrder[0];
            for (size_t n = 0; n < dihedrals.size(); n++)
            {
                double angle = pref.angles[n][metadataIndex];
                cds::AngleWithMetadata am {angle, angle, metadataIndex};
                setDihedralAngle(graph, data, mutableData, dihedrals[n], am);
            }
        };
        std::function<void(const cds::PermutationShapePreference&)> onPermutation =
            [&](const cds::PermutationShapePreference& pref)
        {
            if (rotamerType != GlycamMetadata::RotamerType::permutation)
            {
                throw std::runtime_error(
                    "expected but did not receive PermutationShapePreference for permutation linkage");
            }
            for (size_t n = 0; n < dihedrals.size(); n++)
            {
                size_t metadataIndex = pref.metadataOrder[n][0];
                double angle         = pref.angles[n][metadataIndex];
                cds::AngleWithMetadata am {angle, angle, metadataIndex};
                setDihedralAngle(graph, data, mutableData, dihedrals[n], am);
            }
        };
        return onResidueLinkageShapePreference(onConformer, onPermutation, preference);
    }
} // namespace glycoproteinBuilder
