#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <functional>

namespace glycoproteinBuilder
{
    namespace
    {
        std::vector<bool> groupContains(const std::vector<size_t>& group, const std::vector<size_t>& indices)
        {
            std::vector<bool> result(group.size(), false);
            for (size_t n : indices)
            {
                result[group[n]] = true;
            }
            return result;
        }
    } // namespace

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
        const GlycanIndices& glycan = data.indices.glycans[glycanId];
        size_t siteResidue          = glycan.attachmentResidue;
        updateResidueBounds(graph, mutableData, siteResidue);
        updateResidueMoleculeBounds(graph, mutableData, siteResidue);
        size_t moleculeId = glycan.glycanMolecule;
        for (size_t residue : moleculeResidues(graph, moleculeId))
        {
            updateResidueBounds(graph, mutableData, residue);
        }
        updateMoleculeBounds(graph, mutableData, moleculeId);
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
                          size_t linkageId, size_t dihedralId, const cds::AngleWithMetadata& target)
    {
        const RotatableDihedralIndices& dihedral      = data.indices.rotatableDihedrals[dihedralId];
        std::array<cds::Coordinate, 4> dihedralCoords = dihedralCoordinates(data, mutableData, dihedralId);
        auto matrix = cds::rotationTo(dihedralCoords, constants::toRadians(target.value));
        for (size_t atom : dihedral.movingAtoms)
        {
            cds::Coordinate& coord = mutableData.atomBounds[atom].center;
            coord                  = matrix * coord;
        }
        size_t edgeId                           = data.indices.residueLinkages[linkageId].residueEdge;
        const std::array<size_t, 2>& residueIds = graph.residues.edges.nodeAdjacencies[edgeId];
        std::vector<bool> residueMoving         = groupContains(graph.atomResidue, dihedral.movingAtoms);
        for (size_t n = 0; n < graph.residueCount; n++)
        {
            if (residueMoving[n] && (n != residueIds[0]) && (n != residueIds[1]))
            {
                cds::Coordinate& coord = mutableData.residueBounds[n].center;
                coord                  = matrix * coord;
            }
        }
        for (size_t n : residueIds)
        {
            mutableData.residueBounds[n] =
                cds::boundingSphere(codeUtils::indicesToValues(mutableData.atomBounds, residueAtoms(graph, n)));
        }

        mutableData.currentDihedralShape[dihedralId] = target;
    }

    void setLinkageShape(const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                         size_t glycanId, const std::vector<cds::AngleWithMetadata>& recordedShape)
    {
        const std::vector<size_t>& linkages = data.indices.glycans[glycanId].linkages;
        for (size_t n = 0; n < linkages.size(); n++)
        {
            size_t linkageId                     = linkages[n];
            const std::vector<size_t>& dihedrals = data.indices.residueLinkages[linkageId].rotatableDihedrals;
            for (size_t k = 0; k < dihedrals.size(); k++)
            {
                size_t dihedralId = dihedrals[k];
                setDihedralAngle(graph, data, mutableData, linkageId, dihedralId, recordedShape[dihedralId]);
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
                setDihedralAngle(graph, data, mutableData, linkageId, dihedrals[n], am);
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
                setDihedralAngle(graph, data, mutableData, linkageId, dihedrals[n], am);
            }
        };
        return onResidueLinkageShapePreference(onConformer, onPermutation, preference);
    }
} // namespace glycoproteinBuilder
