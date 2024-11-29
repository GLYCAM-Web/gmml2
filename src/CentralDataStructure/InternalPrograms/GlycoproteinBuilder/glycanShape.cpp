#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
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

    void updateResidueBounds(const AssemblyGraphs& graphs, MutableData& mutableData, size_t index)
    {
        mutableData.residueBounds[index] =
            cds::boundingSphere(codeUtils::indexValues(mutableData.atomBounds, residueAtoms(graphs, index)));
    }

    void updateResidueMoleculeBounds(const AssemblyGraphs& graphs, MutableData& mutableData, size_t index)
    {
        size_t moleculeIndex = graphs.indices.residueMolecule[index];
        mutableData.moleculeBounds[moleculeIndex] =
            cds::boundingSphereIncluding(mutableData.moleculeBounds[moleculeIndex], mutableData.residueBounds[index]);
    }

    void updateMoleculeBounds(const AssemblyGraphs& graphs, MutableData& mutableData, size_t index)
    {
        mutableData.moleculeBounds[index] =
            cds::boundingSphere(codeUtils::indexValues(mutableData.residueBounds, moleculeResidues(graphs, index)));
    }

    void updateGlycanBounds(const AssemblyGraphs& graphs, MutableData& mutableData, size_t glycanId)
    {
        const GlycanIndices& glycan = graphs.indices.glycans[glycanId];
        size_t siteResidue          = glycan.attachmentResidue;
        updateResidueBounds(graphs, mutableData, siteResidue);
        updateResidueMoleculeBounds(graphs, mutableData, siteResidue);
        size_t moleculeId = glycan.glycanMolecule;
        for (size_t residue : moleculeResidues(graphs, moleculeId))
        {
            updateResidueBounds(graphs, mutableData, residue);
        }
        updateMoleculeBounds(graphs, mutableData, moleculeId);
    }

    std::array<cds::Coordinate, 4> dihedralCoordinates(const AssemblyGraphs& graphs, const MutableData& mutableData,
                                                       size_t dihedralId)
    {
        const std::array<size_t, 4>& atoms = graphs.indices.rotatableDihedrals[dihedralId].atoms;
        auto coord                         = [&](size_t n)
        {
            return mutableData.atomBounds[atoms[n]].center;
        };
        return {coord(3), coord(2), coord(1), coord(0)};
    }

    void setDihedralAngle(const AssemblyGraphs& graphs, MutableData& mutableData, size_t linkageId, size_t dihedralId,
                          const cds::AngleWithMetadata& target)
    {
        const RotatableDihedralIndices& dihedral      = graphs.indices.rotatableDihedrals[dihedralId];
        std::array<cds::Coordinate, 4> dihedralCoords = dihedralCoordinates(graphs, mutableData, dihedralId);
        auto matrix = cds::rotationTo(dihedralCoords, constants::toRadians(target.value));
        for (size_t atom : dihedral.movingAtoms)
        {
            cds::Coordinate& coord = mutableData.atomBounds[atom].center;
            coord                  = matrix * coord;
        }
        size_t edgeId                           = graphs.indices.residueLinkages[linkageId].residueEdge;
        const std::array<size_t, 2>& residueIds = graphs.residues.edges.nodeAdjacencies[edgeId];
        std::vector<bool> residueMoving         = groupContains(graphs.indices.atomResidue, dihedral.movingAtoms);
        for (size_t n = 0; n < graphs.indices.residueCount; n++)
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
                cds::boundingSphere(codeUtils::indexValues(mutableData.atomBounds, residueAtoms(graphs, n)));
        }

        mutableData.currentDihedralShape[dihedralId] = target;
    }

    void setLinkageShape(const AssemblyGraphs& graphs, MutableData& mutableData, size_t glycanId,
                         const std::vector<cds::AngleWithMetadata>& recordedShape)
    {
        const std::vector<size_t>& linkages = graphs.indices.glycans[glycanId].linkages;
        for (size_t n = 0; n < linkages.size(); n++)
        {
            size_t linkageId                     = linkages[n];
            const std::vector<size_t>& dihedrals = graphs.indices.residueLinkages[linkageId].rotatableDihedrals;
            for (size_t k = 0; k < dihedrals.size(); k++)
            {
                size_t dihedralId = dihedrals[k];
                setDihedralAngle(graphs, mutableData, linkageId, dihedralId, recordedShape[dihedralId]);
            }
        }
    }

    void setLinkageShapeToPreference(const AssemblyGraphs& graphs, const AssemblyData& data, MutableData& mutableData,
                                     size_t linkageId, const cds::ResidueLinkageShapePreference& preference)
    {
        GlycamMetadata::RotamerType rotamerType = data.residueLinkageData.rotamerTypes[linkageId];
        const std::vector<size_t>& dihedrals    = graphs.indices.residueLinkages[linkageId].rotatableDihedrals;

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
                setDihedralAngle(graphs, mutableData, linkageId, dihedrals[n], am);
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
                setDihedralAngle(graphs, mutableData, linkageId, dihedrals[n], am);
            }
        };
        return onResidueLinkageShapePreference(onConformer, onPermutation, preference);
    }
} // namespace glycoproteinBuilder
