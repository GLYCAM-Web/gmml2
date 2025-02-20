#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <functional>

namespace glycoproteinBuilder
{
    void updateResidueBounds(const assembly::Graph& graph, MutableData& mutableData, size_t index)
    {
        mutableData.residueBounds[index] =
            cds::boundingSphere(codeUtils::indicesToValues(mutableData.atomBounds, residueAtoms(graph, index)));
    }

    void updateResidueMoleculeBounds(const assembly::Graph& graph, MutableData& mutableData, size_t index)
    {
        size_t moleculeIndex = graph.residueMolecule[index];
        mutableData.moleculeBounds[moleculeIndex] =
            cds::boundingSphereIncluding(mutableData.moleculeBounds[moleculeIndex], mutableData.residueBounds[index]);
    }

    void updateMoleculeBounds(const assembly::Graph& graph, MutableData& mutableData, size_t index)
    {
        mutableData.moleculeBounds[index] =
            cds::boundingSphere(codeUtils::indicesToValues(mutableData.residueBounds, moleculeResidues(graph, index)));
    }

    void updateGlycanBounds(const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                            size_t glycanId)
    {
        size_t siteResidue = data.glycans.attachmentResidue[glycanId];
        updateResidueBounds(graph, mutableData, siteResidue);
        updateResidueMoleculeBounds(graph, mutableData, siteResidue);
        size_t moleculeId = data.glycans.moleculeId[glycanId];
        for (size_t residue : moleculeResidues(graph, moleculeId))
        {
            updateResidueBounds(graph, mutableData, residue);
        }
        updateMoleculeBounds(graph, mutableData, moleculeId);
    }

    void updateBoundsContainingAtoms(const assembly::Graph& graph, MutableData& mutableData,
                                     const std::vector<size_t>& atoms)
    {
        std::vector<bool> updateResidue(graph.residueCount, false);
        std::vector<bool> updateMolecule(graph.moleculeCount, false);
        for (size_t atom : atoms)
        {
            updateResidue[graph.atomResidue[atom]] = true;
        }
        for (size_t n : codeUtils::boolsToIndices(updateResidue))
        {
            updateMolecule[graph.residueMolecule[n]] = true;
            updateResidueBounds(graph, mutableData, n);
        }
        for (size_t n : codeUtils::boolsToIndices(updateMolecule))
        {
            updateMoleculeBounds(graph, mutableData, n);
        }
    }

    std::array<cds::Coordinate, 4> dihedralCoordinates(const AssemblyData& data, const MutableData& mutableData,
                                                       size_t dihedralId)
    {
        const std::array<size_t, 4>& atoms = data.indices.rotatableDihedrals[dihedralId].atoms;
        auto coord                         = [&](size_t n)
        {
            return mutableData.atomBounds[atoms[n]].center;
        };
        return {coord(3), coord(2), coord(1), coord(0)};
    }

    void setDihedralAngle(const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                          size_t dihedralId, const cds::AngleWithMetadata& target)
    {
        mutableData.dihedralCurrentMetadata[dihedralId] = target.metadataIndex;
        const RotatableDihedralIndices& dihedral        = data.indices.rotatableDihedrals[dihedralId];
        std::array<cds::Coordinate, 4> dihedralCoords   = dihedralCoordinates(data, mutableData, dihedralId);
        auto matrix = cds::rotationTo(dihedralCoords, constants::toRadians(target.value));
        for (size_t atom : dihedral.movingAtoms)
        {
            cds::Coordinate& coord = mutableData.atomBounds[atom].center;
            coord                  = matrix * coord;
        }
        updateBoundsContainingAtoms(graph, mutableData, dihedral.movingAtoms);
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
