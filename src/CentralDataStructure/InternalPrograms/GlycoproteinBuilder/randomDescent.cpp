#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/randomDescent.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapCount.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/random.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <sstream>

namespace glycoproteinBuilder
{
    namespace
    {
        using GlycamMetadata::DihedralAngleData;
        using GlycamMetadata::DihedralAngleDataVector;
        using GlycamMetadata::RotamerType;

        struct PartialDihedralRotationData
        {
            std::vector<bool> atomMoving;
            std::vector<double> residueWeights;
            std::array<std::vector<size_t>, 2> residueIndices;
        };

        PartialDihedralRotationData toRotationInputData(const AssemblyGraphs& graphs, const AssemblyData& data,
                                                        const MutableData& mutableData, const OverlapWeight& weight,
                                                        size_t glycanId, size_t linkageId, size_t dihedralId)
        {
            const ResidueLinkageIndices& linkage     = graphs.indices.residueLinkages[linkageId];
            const RotatableDihedralIndices& dihedral = graphs.indices.rotatableDihedrals[dihedralId];
            const std::vector<size_t>& movingAtoms   = graphs.indices.rotatableDihedrals[dihedralId].movingAtoms;
            std::vector<bool> atomMoving(graphs.indices.atomCount, false);
            for (size_t atom : movingAtoms)
            {
                atomMoving[atom] = true;
            }
            const std::vector<cds::Sphere>& atomBounds     = mutableData.atomBounds;
            const std::vector<cds::Sphere>& residueBounds  = mutableData.residueBounds;
            const std::vector<cds::Sphere>& moleculeBounds = mutableData.moleculeBounds;
            std::vector<double> residueWeights             = data.residues.overlapWeights;
            cds::Sphere movingAtomBounds = cds::boundingSphere(codeUtils::indicesToValues(atomBounds, movingAtoms));
            Coordinate pointA            = atomBounds[dihedral.atoms[1]].center;
            Coordinate pointB            = atomBounds[dihedral.atoms[2]].center;
            cds::Sphere movementBounds   = boundingSphereCenteredOnLine(movingAtomBounds, pointA, pointB);

            std::vector<size_t> intersectingResidues;
            intersectingResidues.reserve(graphs.indices.residueCount);
            codeUtils::insertInto(intersectingResidues, linkage.nonReducingResidues);
            size_t glycanMolecule = graphs.indices.glycans[glycanId].glycanMolecule;
            for (size_t n : graphs.indices.proteinMolecules)
            {
                if (cds::spheresOverlap(constants::overlapTolerance, movementBounds, moleculeBounds[n]))
                {
                    cds::insertIndicesOfIntersection(intersectingResidues, movementBounds, residueBounds,
                                                     moleculeResidues(graphs, n));
                }
            }
            for (size_t n = 0; n < graphs.indices.glycans.size(); n++)
            {
                size_t otherMolecule = graphs.indices.glycans[n].glycanMolecule;
                if (mutableData.glycanIncluded[n] && (otherMolecule != glycanMolecule) &&
                    cds::spheresOverlap(constants::overlapTolerance, movementBounds, moleculeBounds[otherMolecule]))
                {
                    cds::insertIndicesOfIntersection(intersectingResidues, movementBounds, residueBounds,
                                                     moleculeResidues(graphs, otherMolecule));
                }
            }

            for (size_t n : linkage.reducingResidues)
            {
                residueWeights[n] = 1.0;
            }
            for (size_t n : linkage.nonReducingResidues)
            {
                residueWeights[n] = weight.self;
            }
            return {
                atomMoving, residueWeights, {intersectingResidues, linkage.reducingResidues}
            };
        }

        void wigglePermutationLinkage(const AssemblyGraphs& graphs, const AssemblyData& data, MutableData& mutableData,
                                      size_t glycanId, size_t linkageId, const cds::AngleSearchSettings& settings,
                                      const OverlapWeight& weight,
                                      const cds::PermutationShapePreference& shapePreference)
        {
            const std::vector<size_t>& dihedrals = graphs.indices.residueLinkages[linkageId].rotatableDihedrals;
            const std::vector<cds::DihedralAngleDataVector>& dihedralMetadata = data.rotatableDihedralData.metadata;
            //  Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first
            //  rotatable bond in Asn outwards
            for (size_t rn = 0; rn < dihedrals.size(); rn++)
            {
                size_t n                                         = dihedrals.size() - 1 - rn;
                size_t dihedralId                                = dihedrals[n];
                cds::AngleSearchPreference preference            = {settings.deviation, shapePreference.angles[n],
                                                                    shapePreference.metadataOrder[n]};
                const std::array<cds::Coordinate, 4> coordinates = dihedralCoordinates(graphs, mutableData, dihedralId);
                PartialDihedralRotationData partial =
                    toRotationInputData(graphs, data, mutableData, weight, glycanId, linkageId, dihedralId);
                cds::DihedralRotationData input {partial.atomMoving,
                                                 mutableData.atomBounds,
                                                 mutableData.residueBounds,
                                                 partial.residueWeights,
                                                 graphs.residues.nodes.elements,
                                                 partial.residueIndices,
                                                 data.residueLinkageData.overlapBonds[linkageId]};
                const GlycamMetadata::DihedralAngleDataVector& metadata = dihedralMetadata[dihedralId];
                cds::AngleOverlap best                                  = cds::wiggleUsingRotamers(
                    settings.angles, coordinates, codeUtils::indexVector(metadata), metadata, preference, input);
                setDihedralAngle(graphs, mutableData, linkageId, dihedralId, best.angle);
            }
        }

        void wiggleConformerLinkage(const AssemblyGraphs& graphs, const AssemblyData& data, MutableData& mutableData,
                                    size_t glycanId, size_t linkageId, const cds::AngleSearchSettings& settings,
                                    const OverlapWeight& weight, const cds::ConformerShapePreference& shapePreference)
        {
            const std::vector<size_t>& dihedrals = graphs.indices.residueLinkages[linkageId].rotatableDihedrals;
            const std::vector<cds::DihedralAngleDataVector>& dihedralMetadata = data.rotatableDihedralData.metadata;
            size_t numberOfMetadata                                           = shapePreference.metadataOrder.size();
            const std::vector<std::vector<double>>& preferenceAngles          = shapePreference.angles;
            const std::vector<bool>& isFrozen                                 = shapePreference.isFrozen;
            std::vector<std::vector<cds::AngleWithMetadata>> results;
            results.resize(numberOfMetadata);
            std::vector<cds::AngleOverlap> bestOverlaps;
            bestOverlaps.resize(numberOfMetadata);
            std::vector<size_t> index = codeUtils::indexVector(dihedralMetadata[dihedrals[0]]);
            for (size_t k = 0; k < numberOfMetadata; k++)
            {
                results[k].resize(dihedrals.size());
                std::vector<size_t> order = {shapePreference.metadataOrder[k]};
                cds::ConformerShapePreference pref {isFrozen, preferenceAngles, order};
                setLinkageShapeToPreference(graphs, data, mutableData, linkageId, pref);
                //  Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first
                //  rotatable bond in Asn outwards
                for (size_t rn = 0; rn < dihedrals.size(); rn++)
                {
                    size_t n                              = dihedrals.size() - 1 - rn;
                    size_t dihedralId                     = dihedrals[n];
                    cds::AngleSearchPreference preference = {isFrozen[n] ? 0.0 : settings.deviation,
                                                             preferenceAngles[n], order};
                    cds::DihedralCoordinates coordinates = dihedralCoordinates(graphs, mutableData, dihedralId);

                    PartialDihedralRotationData partial =
                        toRotationInputData(graphs, data, mutableData, weight, glycanId, linkageId, dihedralId);
                    cds::DihedralRotationData input {partial.atomMoving,
                                                     mutableData.atomBounds,
                                                     mutableData.residueBounds,
                                                     partial.residueWeights,
                                                     graphs.residues.nodes.elements,
                                                     partial.residueIndices,
                                                     data.residueLinkageData.overlapBonds[linkageId]};
                    cds::AngleOverlap best = cds::wiggleUsingRotamers(settings.angles, coordinates, index,
                                                                      dihedralMetadata[dihedralId], preference, input);
                    results[k][n]          = best.angle;
                    bestOverlaps[k]        = best;
                    setDihedralAngle(graphs, mutableData, linkageId, dihedralId, best.angle);
                }
            }
            size_t bestIndex = cds::bestOverlapResultIndex(bestOverlaps);
            for (size_t n = 0; n < dihedrals.size(); n++)
            {
                size_t dihedralId                 = dihedrals[n];
                cds::AngleWithMetadata& bestShape = results[bestIndex][n];
                setDihedralAngle(graphs, mutableData, linkageId, dihedralId, bestShape);
            }
        }
    } // namespace

    void wiggleLinkage(const AssemblyGraphs& graphs, const AssemblyData& data, MutableData& mutableData,
                       size_t glycanId, size_t linkageId, const cds::AngleSearchSettings& searchSettings,
                       const OverlapWeight& weight, const cds::ResidueLinkageShapePreference& shapePreference)
    {
        switch (data.residueLinkageData.rotamerTypes[linkageId])
        {
            case RotamerType::permutation:
                {
                    cds::PermutationShapePreference preference =
                        std::get<cds::PermutationShapePreference>(shapePreference);
                    return wigglePermutationLinkage(graphs, data, mutableData, glycanId, linkageId, searchSettings,
                                                    weight, preference);
                }
            case RotamerType::conformer:
                {
                    cds::ConformerShapePreference preference = std::get<cds::ConformerShapePreference>(shapePreference);
                    return wiggleConformerLinkage(graphs, data, mutableData, glycanId, linkageId, searchSettings,
                                                  weight, preference);
                }
        }
        throw std::runtime_error("unhandled linkage shape preference in glycoproteinOverlapResolution wiggleLinkage");
    }

    void wiggleGlycan(const AssemblyGraphs& graphs, const AssemblyData& data, MutableData& mutableData, size_t glycanId,
                      const cds::AngleSearchSettings& searchSettings, const OverlapWeight& weight,
                      const std::vector<cds::ResidueLinkageShapePreference>& preferences)
    {
        const std::vector<size_t>& linkages = graphs.indices.glycans[glycanId].linkages;
        // wiggling twice gives the first linkages a second chance to resolve in a better structure
        for (size_t k = 0; k < 2; k++)
        {
            for (size_t n = 0; n < linkages.size(); n++)
            {
                size_t linkageId = graphs.indices.glycans[glycanId].linkages[n];
                wiggleLinkage(graphs, data, mutableData, glycanId, linkageId, searchSettings, weight, preferences[n]);
            }
        }
        updateGlycanBounds(graphs, mutableData, glycanId);
    }

    GlycoproteinState randomDescent(pcg32& rng, GlycanShapeRandomizer randomizeShape,
                                    const cds::AngleSearchSettings& searchSettings, uint persistCycles,
                                    const OverlapWeight& overlapWeight, const AssemblyGraphs& graphs,
                                    const AssemblyData& data, MutableData& mutableData,
                                    const GlycoproteinState& initialState)
    {
        std::stringstream logss;
        logss << "Random Decent, persisting for " << persistCycles << " cycles.\n";
        gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());

        cds::Overlap globalOverlap       = initialState.overlap;
        std::vector<size_t> overlapSites = initialState.overlapSites;
        auto glycositePreferences        = initialState.preferences;
        uint cycle                       = 0;
        while ((!overlapSites.empty()) && (cycle < persistCycles))
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles));
            cycle++;
            cds::Overlap newGlobalOverlap = globalOverlap;
            for (auto& glycanId : codeUtils::shuffleVector(rng, overlapSites))
            {
                const std::vector<size_t>& linkageIds = graphs.indices.glycans[glycanId].linkages;
                cds::Overlap previousOverlap = localOverlap(graphs, data, mutableData, glycanId, overlapWeight.self);
                auto preferences             = randomizeShape(rng, graphs, data, mutableData, glycanId);
                std::vector<cds::AngleWithMetadata> lastShape = mutableData.currentDihedralShape;
                for (size_t n = 0; n < linkageIds.size(); n++)
                {
                    setLinkageShapeToPreference(graphs, data, mutableData, linkageIds[n], preferences[n]);
                }
                wiggleGlycan(graphs, data, mutableData, glycanId, searchSettings, overlapWeight, preferences);
                updateGlycanBounds(graphs, mutableData, glycanId);
                cds::Overlap newOverlap = localOverlap(graphs, data, mutableData, glycanId, overlapWeight.self);
                cds::Overlap diff       = newOverlap + (previousOverlap * -1);
                bool isWorse            = cds::compareOverlaps(newOverlap, previousOverlap) > 0;
                if (isWorse)
                {
                    setLinkageShape(graphs, mutableData, glycanId, lastShape);
                    updateGlycanBounds(graphs, mutableData, glycanId);
                }
                else
                {
                    newGlobalOverlap               += diff;
                    glycositePreferences[glycanId] = preferences;
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              "RandomDescent accepted a change of " + std::to_string(diff.count));
                }
            }
            if (cds::compareOverlaps(globalOverlap, newGlobalOverlap) > 0)
            {
                cycle = 0;
            }
            globalOverlap = newGlobalOverlap;
            overlapSites  = determineSitesWithOverlap(overlapSites, graphs, data, mutableData);
        }
        return {globalOverlap, overlapSites, glycositePreferences};
    }

} // namespace glycoproteinBuilder
