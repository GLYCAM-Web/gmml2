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
#include "includes/External_Libraries/PCG/pcg_extras.h"

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

        PartialDihedralRotationData toRotationInputData(const AssemblyGraphs& graphs, AssemblyData& data,
                                                        const OverlapWeight& weight, size_t glycanId, size_t linkageId,
                                                        size_t dihedralId)
        {
            const ResidueLinkageIndices& linkage     = graphs.residueLinkages[linkageId];
            const RotatableDihedralIndices& dihedral = graphs.rotatableDihedralIndices[dihedralId];
            const std::vector<size_t>& movingAtoms   = graphs.rotatableDihedralIndices[dihedralId].movingAtoms;
            std::vector<bool> atomMoving(graphs.indices.atoms.size(), false);
            for (size_t atom : movingAtoms)
            {
                atomMoving[atom] = true;
            }
            std::vector<cds::Sphere>& atomBounds    = data.atoms.bounds;
            std::vector<cds::Sphere>& residueBounds = data.residues.bounds;
            std::vector<double> residueWeights      = data.residues.overlapWeights;
            cds::Sphere movingAtomBounds = cds::boundingSphere(codeUtils::indexValues(atomBounds, movingAtoms));
            Coordinate pointA            = data.atoms.bounds[dihedral.atoms[1]].center;
            Coordinate pointB            = data.atoms.bounds[dihedral.atoms[2]].center;
            cds::Sphere movementBounds   = boundingSphereCenteredOnLine(movingAtomBounds, pointA, pointB);

            std::vector<size_t> intersectingResidues;
            intersectingResidues.reserve(graphs.indices.residues.size());
            codeUtils::insertInto(intersectingResidues, linkage.nonReducingResidues);
            size_t glycanMolecule = graphs.glycans[glycanId].glycanMolecule;
            for (size_t n : graphs.proteinMolecules)
            {
                if (cds::spheresOverlap(constants::overlapTolerance, movementBounds, data.molecules.bounds[n]))
                {
                    cds::insertIndicesOfIntersection(intersectingResidues, movementBounds, residueBounds,
                                                     moleculeResidues(graphs, n));
                }
            }
            for (size_t n = 0; n < graphs.glycans.size(); n++)
            {
                size_t otherMolecule = graphs.glycans[n].glycanMolecule;
                if (data.glycanData.included[n] && (otherMolecule != glycanMolecule) &&
                    cds::spheresOverlap(constants::overlapTolerance, movementBounds,
                                        data.molecules.bounds[otherMolecule]))
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

        void wigglePermutationLinkage(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId,
                                      size_t linkageId, const cds::AngleSearchSettings& settings,
                                      const OverlapWeight& weight,
                                      const cds::PermutationShapePreference& shapePreference)
        {
            const std::vector<size_t>& dihedrals       = graphs.residueLinkages[linkageId].rotatableDihedrals;
            const cds::DihedralAngleMetadata& metadata = data.residueLinkageData.metadata[linkageId];
            //  Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first
            //  rotatable bond in Asn outwards
            for (size_t rn = 0; rn < dihedrals.size(); rn++)
            {
                size_t n                                         = dihedrals.size() - 1 - rn;
                size_t dihedralId                                = dihedrals[n];
                cds::AngleSearchPreference preference            = {settings.deviation, shapePreference.angles[n],
                                                                    shapePreference.metadataOrder[n]};
                const std::array<cds::Coordinate, 4> coordinates = dihedralCoordinates(graphs, data, dihedralId);
                PartialDihedralRotationData partial =
                    toRotationInputData(graphs, data, weight, glycanId, linkageId, dihedralId);
                cds::DihedralRotationData input {partial.atomMoving,
                                                 data.atoms.bounds,
                                                 data.residues.bounds,
                                                 partial.residueWeights,
                                                 graphs.residues.nodes.elements,
                                                 partial.residueIndices,
                                                 data.residueLinkageData.overlapBonds[linkageId]};
                const GlycamMetadata::DihedralAngleDataVector& metadataVector = metadata[n];
                cds::AngleOverlap best =
                    cds::wiggleUsingRotamers(settings.angles, coordinates, codeUtils::indexVector(metadataVector),
                                             metadataVector, preference, input);
                setDihedralAngle(graphs, data, linkageId, dihedralId, best.angle);
            }
        }

        void wiggleConformerLinkage(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId, size_t linkageId,
                                    const cds::AngleSearchSettings& settings, const OverlapWeight& weight,
                                    const cds::ConformerShapePreference& shapePreference)
        {
            const std::vector<size_t>& dihedrals       = graphs.residueLinkages[linkageId].rotatableDihedrals;
            const cds::DihedralAngleMetadata& metadata = data.residueLinkageData.metadata[linkageId];
            size_t numberOfMetadata                    = shapePreference.metadataOrder.size();
            const std::vector<std::vector<double>>& preferenceAngles = shapePreference.angles;
            const std::vector<bool>& isFrozen                        = shapePreference.isFrozen;
            std::vector<std::vector<cds::AngleWithMetadata>> results;
            results.resize(numberOfMetadata);
            std::vector<cds::AngleOverlap> bestOverlaps;
            bestOverlaps.resize(numberOfMetadata);
            std::vector<size_t> index = codeUtils::indexVector(metadata[0]);
            for (size_t k = 0; k < numberOfMetadata; k++)
            {
                results[k].resize(dihedrals.size());
                std::vector<size_t> order = {shapePreference.metadataOrder[k]};
                cds::ConformerShapePreference pref {isFrozen, preferenceAngles, order};
                setLinkageShapeToPreference(graphs, data, linkageId, pref);
                //  Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first
                //  rotatable bond in Asn outwards
                for (size_t rn = 0; rn < dihedrals.size(); rn++)
                {
                    size_t n                              = dihedrals.size() - 1 - rn;
                    size_t dihedralId                     = dihedrals[n];
                    cds::AngleSearchPreference preference = {isFrozen[n] ? 0.0 : settings.deviation,
                                                             preferenceAngles[n], order};
                    cds::DihedralCoordinates coordinates = dihedralCoordinates(graphs, data, dihedralId);

                    PartialDihedralRotationData partial =
                        toRotationInputData(graphs, data, weight, glycanId, linkageId, dihedralId);
                    cds::DihedralRotationData input {partial.atomMoving,
                                                     data.atoms.bounds,
                                                     data.residues.bounds,
                                                     partial.residueWeights,
                                                     graphs.residues.nodes.elements,
                                                     partial.residueIndices,
                                                     data.residueLinkageData.overlapBonds[linkageId]};
                    cds::AngleOverlap best =
                        cds::wiggleUsingRotamers(settings.angles, coordinates, index, metadata[n], preference, input);
                    results[k][n]   = best.angle;
                    bestOverlaps[k] = best;
                    setDihedralAngle(graphs, data, linkageId, dihedralId, best.angle);
                }
            }
            size_t bestIndex = cds::bestOverlapResultIndex(bestOverlaps);
            for (size_t n = 0; n < dihedrals.size(); n++)
            {
                size_t dihedralId                 = dihedrals[n];
                cds::AngleWithMetadata& bestShape = results[bestIndex][n];
                setDihedralAngle(graphs, data, linkageId, dihedralId, bestShape);
            }
        }
    } // namespace

    void wiggleLinkage(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId, size_t linkageId,
                       const cds::AngleSearchSettings& searchSettings, const OverlapWeight& weight,
                       const cds::ResidueLinkageShapePreference& shapePreference)
    {
        switch (data.residueLinkageData.rotamerTypes[linkageId])
        {
            case RotamerType::permutation:
                {
                    cds::PermutationShapePreference preference =
                        std::get<cds::PermutationShapePreference>(shapePreference);
                    return wigglePermutationLinkage(graphs, data, glycanId, linkageId, searchSettings, weight,
                                                    preference);
                }
            case RotamerType::conformer:
                {
                    cds::ConformerShapePreference preference = std::get<cds::ConformerShapePreference>(shapePreference);
                    return wiggleConformerLinkage(graphs, data, glycanId, linkageId, searchSettings, weight,
                                                  preference);
                }
        }
        throw std::runtime_error("unhandled linkage shape preference in glycoproteinOverlapResolution wiggleLinkage");
    }

    void wiggleGlycan(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId,
                      const cds::AngleSearchSettings& searchSettings, const OverlapWeight& weight,
                      const std::vector<cds::ResidueLinkageShapePreference>& preferences)
    {
        const std::vector<size_t>& linkages = graphs.glycans[glycanId].linkages;
        // wiggling twice gives the first linkages a second chance to resolve in a better structure
        for (size_t k = 0; k < 2; k++)
        {
            for (size_t n = 0; n < linkages.size(); n++)
            {
                size_t linkageId = graphs.glycans[glycanId].linkages[n];
                wiggleLinkage(graphs, data, glycanId, linkageId, searchSettings, weight, preferences[n]);
            }
        }
        updateGlycanBounds(graphs, data, glycanId);
    }

    GlycoproteinState randomDescent(pcg32 rng, GlycanShapeRandomizer randomizeShape,
                                    const cds::AngleSearchSettings& searchSettings, uint persistCycles,
                                    const OverlapWeight& overlapWeight, const AssemblyGraphs& graphs,
                                    AssemblyData& data, const GlycoproteinState& initialState)
    {
        std::stringstream logss;
        logss << "Random Decent, persisting for " << persistCycles << " cycles.\n";
        gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());

        cds::Overlap globalOverlap                    = initialState.overlap;
        std::vector<size_t> overlapSites              = initialState.overlapSites;
        auto glycositePreferences                     = initialState.preferences;
        std::vector<cds::AngleWithMetadata> lastShape = data.rotatableDihedralData.currentShape;
        uint cycle                                    = 0;
        while ((!overlapSites.empty()) && (cycle < persistCycles))
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles));
            cycle++;
            cds::Overlap newGlobalOverlap = globalOverlap;
            for (auto& glycanId : codeUtils::shuffleVector(rng, overlapSites))
            {
                const std::vector<size_t>& linkageIds = graphs.glycans[glycanId].linkages;
                cds::Overlap previousOverlap          = localOverlap(graphs, data, glycanId, overlapWeight.self);
                auto preferences                      = randomizeShape(graphs, data, glycanId);
                lastShape                             = data.rotatableDihedralData.currentShape;
                for (size_t n = 0; n < linkageIds.size(); n++)
                {
                    setLinkageShapeToPreference(graphs, data, linkageIds[n], preferences[n]);
                }
                wiggleGlycan(graphs, data, glycanId, searchSettings, overlapWeight, preferences);
                updateGlycanBounds(graphs, data, glycanId);
                cds::Overlap newOverlap = localOverlap(graphs, data, glycanId, overlapWeight.self);
                cds::Overlap diff       = newOverlap + (previousOverlap * -1);
                bool isWorse            = cds::compareOverlaps(newOverlap, previousOverlap) > 0;
                if (isWorse)
                {
                    setLinkageShape(graphs, data, glycanId, lastShape);
                    updateGlycanBounds(graphs, data, glycanId);
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
            overlapSites  = determineSitesWithOverlap(overlapSites, graphs, data);
        }
        return {globalOverlap, overlapSites, glycositePreferences};
    }

} // namespace glycoproteinBuilder
