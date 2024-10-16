#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"

#include "includes/CodeUtils/containers.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"

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

    void updateResidueBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t index)
    {
        data.residues.bounds[index] =
            cds::boundingSphere(codeUtils::indexValues(data.atoms.bounds, residueAtoms(graphs, index)));
    }

    void updateResidueMoleculeBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t index)
    {
        size_t moleculeIndex = graphs.indices.residueMolecule[index];
        data.molecules.bounds[moleculeIndex] =
            cds::boundingSphereIncluding(data.molecules.bounds[moleculeIndex], data.residues.bounds[index]);
    }

    void updateMoleculeBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t index)
    {
        data.molecules.bounds[index] =
            cds::boundingSphere(codeUtils::indexValues(data.residues.bounds, moleculeResidues(graphs, index)));
    }

    void updateGlycanBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId)
    {
        const GlycanIndices& glycan = graphs.glycans[glycanId];
        size_t siteResidue          = glycan.attachmentResidue;
        updateResidueBounds(graphs, data, siteResidue);
        updateResidueMoleculeBounds(graphs, data, siteResidue);
        size_t moleculeId = glycan.glycanMolecule;
        for (size_t residue : moleculeResidues(graphs, moleculeId))
        {
            updateResidueBounds(graphs, data, residue);
        }
        updateMoleculeBounds(graphs, data, moleculeId);
    }

    std::array<cds::Coordinate, 4> dihedralCoordinates(const AssemblyGraphs& graphs, const AssemblyData& data,
                                                       size_t dihedralId)
    {
        const std::array<size_t, 4>& atoms = graphs.rotatableDihedralIndices[dihedralId].atoms;
        auto coord                         = [&](size_t n)
        {
            return data.atoms.bounds[atoms[n]].center;
        };
        return {coord(3), coord(2), coord(1), coord(0)};
    }

    void setDihedralAngle(const AssemblyGraphs& graphs, AssemblyData& data, size_t linkageId, size_t dihedralId,
                          const cds::AngleWithMetadata& target)
    {
        const RotatableDihedralIndices& dihedral      = graphs.rotatableDihedralIndices[dihedralId];
        std::array<cds::Coordinate, 4> dihedralCoords = dihedralCoordinates(graphs, data, dihedralId);
        auto matrix = cds::rotationTo(dihedralCoords, constants::toRadians(target.value));
        for (size_t atom : dihedral.movingAtoms)
        {
            cds::Coordinate& coord = data.atoms.bounds[atom].center;
            coord                  = matrix * coord;
        }
        size_t edgeId                           = graphs.residueLinkages[linkageId].residueEdge;
        const std::array<size_t, 2>& residueIds = graphs.residues.edges.nodeAdjacencies[edgeId];
        std::vector<bool> residueMoving         = groupContains(graphs.indices.atomResidue, dihedral.movingAtoms);
        for (size_t n = 0; n < graphs.indices.residues.size(); n++)
        {
            if (residueMoving[n] && (n != residueIds[0]) && (n != residueIds[1]))
            {
                cds::Coordinate& coord = data.residues.bounds[n].center;
                coord                  = matrix * coord;
            }
        }
        for (size_t n : residueIds)
        {
            data.residues.bounds[n] =
                cds::boundingSphere(codeUtils::indexValues(data.atoms.bounds, residueAtoms(graphs, n)));
        }

        data.rotatableDihedralData.currentShape[dihedralId] = target;
    }

    void setLinkageShape(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId,
                         const std::vector<cds::AngleWithMetadata>& recordedShape)
    {
        const std::vector<size_t>& linkages = graphs.glycans[glycanId].linkages;
        for (size_t n = 0; n < linkages.size(); n++)
        {
            size_t linkageId                     = linkages[n];
            const std::vector<size_t>& dihedrals = graphs.residueLinkages[linkageId].rotatableDihedrals;
            for (size_t k = 0; k < dihedrals.size(); k++)
            {
                size_t dihedralId = dihedrals[k];
                setDihedralAngle(graphs, data, linkageId, dihedralId, recordedShape[dihedralId]);
            }
        }
    }

    void setLinkageShapeToPreference(const AssemblyGraphs& graphs, AssemblyData& data, size_t linkageId,
                                     const cds::ResidueLinkageShapePreference& preference)
    {
        const std::vector<size_t>& dihedrals = graphs.residueLinkages[linkageId].rotatableDihedrals;
        if (data.residueLinkageData.rotamerTypes[linkageId] == GlycamMetadata::RotamerType::conformer)
        {
            if (!std::holds_alternative<cds::ConformerShapePreference>(preference))
            {
                throw std::runtime_error("expected but did not receive ConformerShapePreference for conformer linkage");
            }
            cds::ConformerShapePreference pref = std::get<cds::ConformerShapePreference>(preference);
            size_t metadataIndex               = pref.metadataOrder[0];
            for (size_t n = 0; n < dihedrals.size(); n++)
            {
                double angle = pref.angles[n][metadataIndex];
                cds::AngleWithMetadata am {angle, angle, metadataIndex};
                setDihedralAngle(graphs, data, linkageId, dihedrals[n], am);
            }
        }
        else
        {
            if (!std::holds_alternative<cds::PermutationShapePreference>(preference))
            {
                throw std::runtime_error(
                    "expected but did not receive PermutationShapePreference for permutation linkage");
            }
            cds::PermutationShapePreference pref = std::get<cds::PermutationShapePreference>(preference);
            for (size_t n = 0; n < dihedrals.size(); n++)
            {
                size_t metadataIndex = pref.metadataOrder[n][0];
                double angle         = pref.angles[n][metadataIndex];
                cds::AngleWithMetadata am {angle, angle, metadataIndex};
                setDihedralAngle(graphs, data, linkageId, dihedrals[n], am);
            }
        }
    }
} // namespace glycoproteinBuilder
