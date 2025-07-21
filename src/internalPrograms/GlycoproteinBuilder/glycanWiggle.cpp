#include "include/internalPrograms/GlycoproteinBuilder/glycanWiggle.hpp"

#include "include/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "include/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblySelection.hpp"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "include/util/containers.hpp"

#include <sstream>

namespace gmml
{
    namespace gpbuilder
    {
        namespace
        {
            void wigglePermutationLinkage(
                const assembly::Graph& graph,
                const AssemblyData& data,
                const assembly::Selection& selection,
                MutableData& mutableData,
                size_t linkageId,
                const AngleSearchSettings& settings,
                const PermutationShapePreference& shapePreference)
            {
                const std::vector<size_t>& dihedrals = data.indices.residueLinkages[linkageId].rotatableDihedrals;
                const std::vector<std::vector<size_t>>& dihedralMetadata = data.rotatableDihedralData.metadata;
                //  Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first
                //  rotatable bond in Asn outwards
                for (size_t rn = 0; rn < dihedrals.size(); rn++)
                {
                    size_t n = dihedrals.size() - 1 - rn;
                    size_t dihedralId = dihedrals[n];
                    AngleSearchPreference preference = {
                        settings.deviation, shapePreference.angles[n], shapePreference.metadataOrder[n]};
                    const std::array<Coordinate, 4> coordinates =
                        dihedralCoordinates(data, mutableData.bounds, dihedralId);
                    const std::vector<size_t>& movingAtoms = data.indices.rotatableDihedrals[dihedralId].movingAtoms;
                    std::vector<bool> atomMoving = util::indicesToBools(graph.indices.atomCount, movingAtoms);
                    assembly::Selection moving =
                        assembly::intersection(graph, selection, assembly::selectByAtoms(graph, atomMoving));
                    assembly::Selection nonMoving = assembly::intersection(
                        graph, selection, assembly::selectByAtoms(graph, util::vectorNot(atomMoving)));

                    const std::vector<size_t>& metadata = dihedralMetadata[dihedralId];
                    auto searchOverlap = [&](const assembly::Bounds& bounds)
                    {
                        return overlapAboveThresholdSum(
                            data.overlapRejectionThreshold,
                            overlapsBetweenSelections(
                                data.potentialTable,
                                data.overlapTolerance,
                                graph,
                                bounds,
                                moving,
                                nonMoving,
                                data.atoms.elements,
                                data.residueEdges.atomsCloseToEdge));
                    };
                    OverlapState best = wiggleUsingRotamers(
                        searchOverlap,
                        settings.angles,
                        data.dihedralAngleTable,
                        graph,
                        mutableData.bounds,
                        movingAtoms,
                        coordinates,
                        util::indexVector(metadata),
                        metadata,
                        preference);
                    mutableData.bounds = best.bounds;
                    mutableData.dihedralCurrentMetadata[dihedralId] = best.angle.metadataIndex;
                }
            }

            void wiggleConformerLinkage(
                const assembly::Graph& graph,
                const AssemblyData& data,
                const assembly::Selection& selection,
                MutableData& mutableData,
                size_t linkageId,
                const AngleSearchSettings& settings,
                const ConformerShapePreference& shapePreference)
            {
                const std::vector<size_t>& dihedrals = data.indices.residueLinkages[linkageId].rotatableDihedrals;
                const std::vector<std::vector<size_t>>& dihedralMetadata = data.rotatableDihedralData.metadata;
                size_t numberOfMetadata = shapePreference.metadataOrder.size();
                const std::vector<std::vector<double>>& preferenceAngles = shapePreference.angles;
                const std::vector<bool>& isFrozen = shapePreference.isFrozen;
                std::vector<assembly::Selection> dihedralMoving;
                dihedralMoving.reserve(dihedrals.size());
                std::vector<assembly::Selection> dihedralNonMoving;
                dihedralNonMoving.reserve(dihedrals.size());
                for (size_t rn = 0; rn < dihedrals.size(); rn++)
                {
                    size_t n = dihedrals.size() - 1 - rn;
                    size_t dihedralId = dihedrals[n];
                    const std::vector<size_t>& movingAtoms = data.indices.rotatableDihedrals[dihedralId].movingAtoms;
                    std::vector<bool> atomMoving = util::indicesToBools(graph.indices.atomCount, movingAtoms);
                    assembly::Selection moving =
                        assembly::intersection(graph, selection, assembly::selectByAtoms(graph, atomMoving));
                    assembly::Selection nonMoving = assembly::intersection(
                        graph, selection, assembly::selectByAtoms(graph, util::vectorNot(atomMoving)));
                    dihedralMoving.push_back(moving);
                    dihedralNonMoving.push_back(nonMoving);
                }

                std::vector<AngleOverlap> bestResults;
                bestResults.resize(numberOfMetadata + 1);
                std::vector<assembly::Bounds> states;
                states.resize(numberOfMetadata + 1);
                std::vector<size_t> index = util::indexVector(dihedralMetadata[dihedrals[0]]);
                assembly::Bounds initialBounds = mutableData.bounds;
                auto runIterationWithMetadata = [&](size_t iteration, size_t metadataIndex)
                {
                    std::vector<size_t> order = {metadataIndex};
                    ConformerShapePreference pref {isFrozen, preferenceAngles, order};
                    //  Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first
                    //  rotatable bond in Asn outwards
                    for (size_t rn = 0; rn < dihedrals.size(); rn++)
                    {
                        size_t n = dihedrals.size() - 1 - rn;
                        size_t dihedralId = dihedrals[n];
                        AngleSearchPreference preference = {
                            isFrozen[n] ? 0.0 : settings.deviation, preferenceAngles[n], order};
                        DihedralCoordinates coordinates = dihedralCoordinates(data, mutableData.bounds, dihedralId);
                        const std::vector<size_t>& movingAtoms =
                            data.indices.rotatableDihedrals[dihedralId].movingAtoms;

                        auto searchOverlap = [&](const assembly::Bounds& bounds)
                        {
                            return overlapAboveThresholdSum(
                                data.overlapRejectionThreshold,
                                overlapsBetweenSelections(
                                    data.potentialTable,
                                    data.overlapTolerance,
                                    graph,
                                    bounds,
                                    dihedralMoving[n],
                                    dihedralNonMoving[n],
                                    data.atoms.elements,
                                    data.residueEdges.atomsCloseToEdge));
                        };
                        OverlapState best = wiggleUsingRotamers(
                            searchOverlap,
                            settings.angles,
                            data.dihedralAngleTable,
                            graph,
                            mutableData.bounds,
                            movingAtoms,
                            coordinates,
                            index,
                            dihedralMetadata[dihedralId],
                            preference);
                        mutableData.bounds = best.bounds;
                    }
                    std::vector<double> overlap = overlapsBetweenSelections(
                        data.potentialTable,
                        data.overlapTolerance,
                        graph,
                        mutableData.bounds,
                        dihedralMoving[0],
                        dihedralNonMoving[0],
                        data.atoms.elements,
                        data.residueEdges.atomsCloseToEdge);
                    bestResults[iteration] = {
                        overlapVectorSum(overlap), {0.0, 0.0, metadataIndex}
                    };
                    states[iteration] = mutableData.bounds;
                };
                size_t initialMetadata = mutableData.dihedralCurrentMetadata[dihedrals[0]];
                // start by looking at current shape with current metadata
                runIterationWithMetadata(0, initialMetadata);
                // then look at the default shape for each metadata
                for (size_t k = 0; k < numberOfMetadata; k++)
                {
                    size_t usedMetadata = shapePreference.metadataOrder[k];
                    std::vector<size_t> order = {usedMetadata};
                    ConformerShapePreference pref {isFrozen, preferenceAngles, order};
                    setLinkageShapeToPreference(graph, data, mutableData, linkageId, pref);
                    runIterationWithMetadata(k + 1, usedMetadata);
                }
                size_t bestIndex = bestOverlapResultIndex(bestResults);
                size_t bestMetadata = bestResults[bestIndex].angle.metadataIndex;
                mutableData.bounds = states[bestIndex];
                for (size_t n = 0; n < dihedrals.size(); n++)
                {
                    mutableData.dihedralCurrentMetadata[dihedrals[n]] = bestMetadata;
                }
            }
        } // namespace

        void wiggleLinkage(
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            MutableData& mutableData,
            size_t linkageId,
            const AngleSearchSettings& searchSettings,
            const ResidueLinkageShapePreference& shapePreference)
        {
            switch (data.residueLinkageData.rotamerTypes[linkageId])
            {
                case RotamerType::permutation:
                    {
                        PermutationShapePreference preference = std::get<PermutationShapePreference>(shapePreference);
                        return wigglePermutationLinkage(
                            graph, data, selection, mutableData, linkageId, searchSettings, preference);
                    }
                case RotamerType::conformer:
                    {
                        ConformerShapePreference preference = std::get<ConformerShapePreference>(shapePreference);
                        return wiggleConformerLinkage(
                            graph, data, selection, mutableData, linkageId, searchSettings, preference);
                    }
            }
            throw std::runtime_error(
                "unhandled linkage shape preference in glycoproteinOverlapResolution wiggleLinkage");
        }

        void wiggleGlycan(
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            const AngleSearchSettings& searchSettings,
            const GlycanShapePreference& preferences,
            MutableData& mutableData,
            size_t glycanId)
        {
            const std::vector<size_t>& linkages = data.glycans.linkages[glycanId];
            // wiggling twice gives the first linkages a second chance to resolve in a better structure
            for (size_t k = 0; k < 2; k++)
            {
                for (size_t n = 0; n < linkages.size(); n++)
                {
                    size_t linkageId = linkages[n];
                    wiggleLinkage(graph, data, selection, mutableData, linkageId, searchSettings, preferences[n]);
                }
            }
        }
    } // namespace gpbuilder
} // namespace gmml
