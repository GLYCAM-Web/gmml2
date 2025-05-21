#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/randomDescent.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapCount.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanWiggle.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/Assembly/assemblySelection.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/random.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <sstream>

namespace glycoproteinBuilder
{
    GlycoproteinState randomDescent(pcg32& rng, const AngleSettings& settings, WiggleGlycan wiggleGlycan,
                                    GlycanShapeRandomizer randomizeShape, SidechainAdjustment adjustSidechains,
                                    uint persistCycles, const assembly::Graph& graph, const AssemblyData& data,
                                    MutableData& mutableData, const GlycoproteinState& initialState)
    {
        std::stringstream logss;
        logss << "Random Decent, persisting for " << persistCycles << " cycles.\n";
        gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());

        cds::Overlap globalOverlap                                   = initialState.overlap;
        std::vector<size_t> sitesWithOverlap                         = initialState.sitesWithOverlap;
        std::vector<size_t> sitesAboveOverlapThreshold               = initialState.sitesAboveOverlapThreshold;
        std::vector<cds::GlycanShapePreference> glycositePreferences = initialState.preferences;
        const assembly::Selection eachSelection                      = assembly::selectByAtomsAndMolecules(
            graph, data.atoms.includeInEachOverlapCheck, mutableData.moleculeIncluded);
        const assembly::Selection mainSelection = assembly::selectByAtomsAndMolecules(
            graph, data.atoms.includeInMainOverlapCheck, mutableData.moleculeIncluded);
        uint cycle = 0;
        while ((!sitesAboveOverlapThreshold.empty()) && (cycle < persistCycles))
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles));
            cycle++;
            cds::Overlap newGlobalOverlap = globalOverlap;
            for (auto& glycanId : codeUtils::shuffleVector(rng, sitesWithOverlap))
            {
                const std::vector<size_t>& linkageIds = data.glycans.linkages[glycanId];
                cds::Overlap previousOverlap = localOverlap(graph, data, eachSelection, mutableData.bounds, glycanId);
                std::vector<cds::GlycanShapePreference> currentPreferences = glycositePreferences;
                currentPreferences[glycanId] = randomizeShape(rng, settings, data, mutableData, glycanId);
                cds::GlycanShapePreference& glycanPreferences = currentPreferences[glycanId];
                MutableData lastShape                         = mutableData;
                for (size_t n = 0; n < linkageIds.size(); n++)
                {
                    setLinkageShapeToPreference(graph, data, mutableData, linkageIds[n], glycanPreferences[n]);
                }
                wiggleGlycan(graph, data, mainSelection, settings, glycanPreferences, mutableData, glycanId);
                adjustSidechains(rng, settings, wiggleGlycan, graph, data, mutableData, currentPreferences, {glycanId});
                cds::Overlap newOverlap = localOverlap(graph, data, eachSelection, mutableData.bounds, glycanId);
                cds::Overlap diff       = newOverlap + (previousOverlap * -1);
                bool isWorse            = cds::compareOverlaps(newOverlap, previousOverlap) > 0;
                if (isWorse)
                {
                    mutableData = lastShape;
                }
                else
                {
                    newGlobalOverlap     += diff;
                    glycositePreferences = currentPreferences;
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              "RandomDescent accepted a change of " + std::to_string(diff.count));
                }
            }

            if (cds::compareOverlaps(globalOverlap, newGlobalOverlap) > 0)
            {
                cycle = 0;
            }
            globalOverlap = newGlobalOverlap;
            sitesWithOverlap =
                determineSitesWithOverlap(0.0, sitesWithOverlap, graph, data, eachSelection, mutableData.bounds);
            sitesAboveOverlapThreshold = determineSitesWithOverlap(data.overlapRejectionThreshold, sitesWithOverlap,
                                                                   graph, data, eachSelection, mutableData.bounds);
        }
        return {globalOverlap, sitesWithOverlap, sitesAboveOverlapThreshold, glycositePreferences};
    }

    void resolveOverlapsWithWiggler(pcg32& rng, const AngleSettings& initialAngleSettings,
                                    const AngleSettings& mainAngleSettings, SidechainAdjustment adjustSidechains,
                                    SidechainAdjustment restoreSidechains, GlycanShapeRandomizer& randomizeShape,
                                    WiggleGlycan wiggleGlycan, const assembly::Graph& graph, const AssemblyData& data,
                                    MutableData& mutableData, size_t persistCycles, bool deleteSitesUntilResolved)
    {
        std::vector<cds::GlycanShapePreference> glycositePreferences;
        const std::vector<size_t> glycanIndices = codeUtils::indexVector(data.glycans.moleculeId);
        for (size_t glycanId : glycanIndices)
        {
            auto preference = randomizeShape(rng, initialAngleSettings, data, mutableData, glycanId);
            const std::vector<size_t>& linkageIds = data.glycans.linkages[glycanId];
            for (size_t k = 0; k < linkageIds.size(); k++)
            {
                setLinkageShapeToPreference(graph, data, mutableData, linkageIds[k], preference[k]);
            }
            glycositePreferences.push_back(preference);
        }
        assembly::Selection glycanSelection = assembly::selectByAtoms(graph, data.atoms.includeInMainOverlapCheck);
        for (size_t glycanId : codeUtils::shuffleVector(rng, glycanIndices))
        {
            wiggleGlycan(graph, data, glycanSelection, initialAngleSettings, glycositePreferences[glycanId],
                         mutableData, glycanId);
        }
        adjustSidechains(rng, initialAngleSettings, wiggleGlycan, graph, data, mutableData, glycositePreferences,
                         glycanIndices);
        assembly::Selection selection = assembly::selectByAtoms(graph, data.atoms.includeInEachOverlapCheck);
        GlycoproteinState currentState;
        std::vector<size_t> sitesWithOverlap =
            determineSitesWithOverlap(0.0, glycanIndices, graph, data, selection, mutableData.bounds);
        std::vector<size_t> sitesAboveOverlapThreshold = determineSitesWithOverlap(
            data.overlapRejectionThreshold, sitesWithOverlap, graph, data, selection, mutableData.bounds);
        for (bool done = false; !done; done = sitesAboveOverlapThreshold.empty() || !deleteSitesUntilResolved)
        {
            cds::Overlap initialOverlap =
                cds::overlapVectorSum(totalOverlaps(graph, data, selection, mutableData.bounds));

            GlycoproteinState initialState = {initialOverlap, sitesWithOverlap, sitesAboveOverlapThreshold,
                                              glycositePreferences};
            currentState     = randomDescent(rng, mainAngleSettings, wiggleGlycan, randomizeShape, adjustSidechains,
                                             persistCycles, graph, data, mutableData, initialState);
            sitesWithOverlap = currentState.sitesWithOverlap;
            sitesAboveOverlapThreshold = currentState.sitesAboveOverlapThreshold;
            if (deleteSitesUntilResolved && !sitesAboveOverlapThreshold.empty())
            {
                size_t indexToRemove = codeUtils::randomIndex(rng, sitesAboveOverlapThreshold);
                size_t glycan        = sitesAboveOverlapThreshold[indexToRemove];
                deleteMolecule(mutableData, data.glycans.moleculeId[glycan]);
                selection             = assembly::selectByAtomsAndMolecules(graph, data.atoms.includeInEachOverlapCheck,
                                                                            mutableData.moleculeIncluded);
                size_t proteinResidue = data.glycans.attachmentResidue[glycan];
                // restore atoms to initial shape
                for (size_t n : residueAtoms(graph, proteinResidue))
                {
                    mutableData.bounds.atoms[n] = data.atoms.initialState[n];
                }
                updateResidueBounds(graph, mutableData.bounds, proteinResidue);
                updateResidueMoleculeBounds(graph, mutableData.bounds, proteinResidue);
                sitesWithOverlap =
                    determineSitesWithOverlap(0.0, includedGlycanIndices(data, mutableData.moleculeIncluded), graph,
                                              data, selection, mutableData.bounds);
                sitesAboveOverlapThreshold = determineSitesWithOverlap(data.overlapRejectionThreshold, sitesWithOverlap,
                                                                       graph, data, selection, mutableData.bounds);
            }
        }
        restoreSidechains(rng, mainAngleSettings, wiggleGlycan, graph, data, mutableData, glycositePreferences,
                          includedGlycanIndices(data, mutableData.moleculeIncluded));
        gmml::log(__LINE__, __FILE__, gmml::INF, "Overlap: " + std::to_string(currentState.overlap.count));
    };
} // namespace glycoproteinBuilder
