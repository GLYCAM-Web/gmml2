#include "include/programs/GlycoproteinBuilder/glycanShape.hpp"

#include "include/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "include/CentralDataStructure/atom.hpp"
#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/util/containers.hpp"

#include <functional>

namespace gmml
{
    namespace gpbuilder
    {
        std::array<Coordinate, 4> dihedralCoordinates(
            const AssemblyData& data, const assembly::Bounds& bounds, size_t dihedralId)
        {
            const std::array<size_t, 4>& atoms = data.indices.rotatableDihedrals[dihedralId].atoms;
            auto coord = [&](size_t n) { return bounds.atoms[atoms[n]].center; };
            return {coord(3), coord(2), coord(1), coord(0)};
        }

        AngleWithMetadata currentDihedralShape(
            const AssemblyData& data, const MutableData& mutableData, size_t dihedralId)
        {
            size_t metadataIndex = mutableData.dihedralCurrentMetadata[dihedralId];
            auto& currentMetadata = data.rotatableDihedralData.metadata[dihedralId][metadataIndex];
            return {
                constants::toDegrees(angle(dihedralCoordinates(data, mutableData.bounds, dihedralId))),
                data.dihedralAngleTable.entries[currentMetadata].default_angle,
                metadataIndex};
        }

        void setDihedralAngle(
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            size_t dihedralId,
            const AngleWithMetadata& target)
        {
            mutableData.dihedralCurrentMetadata[dihedralId] = target.metadataIndex;
            const RotatableDihedralIndices& dihedral = data.indices.rotatableDihedrals[dihedralId];
            std::array<Coordinate, 4> dihedralCoords = dihedralCoordinates(data, mutableData.bounds, dihedralId);
            auto matrix = rotationTo(dihedralCoords, constants::toRadians(target.value));
            for (size_t atom : dihedral.movingAtoms)
            {
                Coordinate& coord = mutableData.bounds.atoms[atom].center;
                coord = matrix * coord;
            }
            updateBoundsContainingAtoms(graph, mutableData.bounds, dihedral.movingAtoms);
        }

        void setLinkageShape(
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            size_t glycanId,
            const std::vector<AngleWithMetadata>& recordedShape)
        {
            const std::vector<size_t>& linkages = data.glycans.linkages[glycanId];
            for (size_t n = 0; n < linkages.size(); n++)
            {
                size_t linkageId = linkages[n];
                const std::vector<size_t>& dihedrals = data.indices.residueLinkages[linkageId].rotatableDihedrals;
                for (size_t k = 0; k < dihedrals.size(); k++)
                {
                    size_t dihedralId = dihedrals[k];
                    setDihedralAngle(graph, data, mutableData, dihedralId, recordedShape[dihedralId]);
                }
            }
        }

        void setLinkageShapeToPreference(
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            size_t linkageId,
            const ResidueLinkageShapePreference& preference)
        {
            RotamerType rotamerType = data.residueLinkageData.rotamerTypes[linkageId];
            const std::vector<size_t>& dihedrals = data.indices.residueLinkages[linkageId].rotatableDihedrals;

            std::function<void(const ConformerShapePreference&)> onConformer = [&](const ConformerShapePreference& pref)
            {
                if (rotamerType != RotamerType::conformer)
                {
                    throw std::runtime_error(
                        "expected but did not receive ConformerShapePreference for conformer linkage");
                }
                size_t metadataIndex = pref.metadataOrder[0];
                for (size_t n = 0; n < dihedrals.size(); n++)
                {
                    double angle = pref.angles[n][metadataIndex];
                    AngleWithMetadata am {angle, angle, metadataIndex};
                    setDihedralAngle(graph, data, mutableData, dihedrals[n], am);
                }
            };
            std::function<void(const PermutationShapePreference&)> onPermutation =
                [&](const PermutationShapePreference& pref)
            {
                if (rotamerType != RotamerType::permutation)
                {
                    throw std::runtime_error(
                        "expected but did not receive PermutationShapePreference for permutation linkage");
                }
                for (size_t n = 0; n < dihedrals.size(); n++)
                {
                    size_t metadataIndex = pref.metadataOrder[n][0];
                    double angle = pref.angles[n][metadataIndex];
                    AngleWithMetadata am {angle, angle, metadataIndex};
                    setDihedralAngle(graph, data, mutableData, dihedrals[n], am);
                }
            };
            return onResidueLinkageShapePreference(onConformer, onPermutation, preference);
        }
    } // namespace gpbuilder
} // namespace gmml
