#include "include/programs/GlycoproteinBuilder/randomDescent.hpp"

#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblySelection.hpp"
#include "include/carbohydrate/dihedralAngleSearch.hpp"
#include "include/carbohydrate/dihedralShape.hpp"
#include "include/external/pcg/pcg_random.h"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/glycoprotein/glycanShape.hpp"
#include "include/glycoprotein/glycanWiggle.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/glycoprotein/glycoproteinUtil.hpp"
#include "include/glycoprotein/overlapCount.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/random.hpp"

#include <sstream>

namespace gmml
{
    namespace gpbuilder
    {
        namespace
        {
            MutableData adjustGlycans(
                pcg32& rng,
                const DihedralAngleDataTable& dihedralAngleTable,
                const OverlapSettings& overlapSettings,
                const AngleSettings& settings,
                size_t iterations,
                const assembly::Graph& graph,
                const assembly::Selection& selection,
                const AssemblyData& data,
                const std::vector<GlycanShapePreference>& shapePreferences,
                const std::vector<size_t>& glycanIds,
                const MutableData& initialState)
            {
                auto standardDeviation = [](const AngleSettings& settings, const DihedralAngleData& metadata)
                {
                    std::function<std::pair<double, double>(const AngleLimit&)> onLimit = [&](const AngleLimit& dev)
                    {
                        double max_std = settings.preferenceDeviation + settings.searchDeviation;
                        double lower_std = dev.lowerDeviationLimit / max_std;
                        double upper_std = dev.upperDeviationLimit / max_std;
                        return std::pair<double, double> {lower_std, upper_std};
                    };
                    std::function<std::pair<double, double>(const AngleStd&)> onStd = [&](const AngleStd& dev) {
                        return std::pair<double, double> {dev.lowerDeviationStd, dev.upperDeviationStd};
                    };
                    return onAngleDeviation(onLimit, onStd, metadata.angle_deviation);
                };
                auto searchAngles = [&standardDeviation, &settings](
                                        const DihedralAngleData& metadata,
                                        double preference,
                                        double deviation,
                                        uint halfIntervalSearches)
                {
                    auto std = standardDeviation(settings, metadata);
                    return AngleSpacing {
                        preference,
                        deviation * std.first,
                        deviation * std.second,
                        settings.searchIncrement,
                        halfIntervalSearches};
                };
                MutableData currentState = initialState;
                size_t maxLinkageCount = 0;
                for (size_t glycanId : glycanIds)
                {
                    maxLinkageCount = std::max(maxLinkageCount, data.glycans.linkages[glycanId].size());
                }
                for (size_t n = 0; n < iterations; n++)
                {
                    for (size_t k = 0; k < maxLinkageCount; k++)
                    {
                        for (size_t glycanId : util::shuffleVector(rng, glycanIds))
                        {
                            const std::vector<size_t>& linkages = data.glycans.linkages[glycanId];
                            if (k < linkages.size())
                            {
                                wiggleLinkage(
                                    dihedralAngleTable,
                                    overlapSettings,
                                    graph,
                                    data,
                                    selection,
                                    currentState,
                                    linkages[k],
                                    {settings.searchDeviation, settings.halfIntervalSearches, searchAngles},
                                    shapePreferences[glycanId][k]);
                            }
                        }
                    }
                }
                return currentState;
            }

            std::vector<std::vector<size_t>> overlapConcerts(const OverlapSites& sites)
            {
                size_t siteCount = sites.concertIds.size();
                std::vector<std::vector<size_t>> result;
                result.reserve(siteCount);
                for (size_t n = 0; n < siteCount; n++)
                {
                    std::vector<size_t> concertGlycans = util::indicesOfElement(sites.concertIds, n);
                    std::vector<bool> hasOverlap = util::indicesToValues(sites.aboveOverlapThreshold, concertGlycans);
                    std::vector<size_t> concert = util::boolsToValues(concertGlycans, hasOverlap);
                    if (concert.size() > 0)
                    {
                        result.push_back(concert);
                    }
                }
                return result;
            }
        } // namespace

        GlycoproteinState randomDescent(
            pcg32& rng,
            const DihedralAngleDataTable& dihedralAngleDataTable,
            PersistCycleAngleSettings toAngleSettings,
            GlycanShapeRandomizer randomizeShape,
            SidechainAdjustment adjustSidechains,
            uint persistCycles,
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const GlycoproteinState& initialState)
        {
            std::stringstream logss;
            logss << "Random Decent, persisting for " << persistCycles << " cycles.\n";
            util::log(__LINE__, __FILE__, util::INF, logss.str());

            const assembly::Selection mainSelection = assembly::selectByAtomsAndMolecules(
                graph, data.atoms.includeInMainOverlapCheck, initialState.mutableData.moleculeIncluded);
            const assembly::Selection fullSelection = assembly::selectByAtomsAndMolecules(
                graph, data.atoms.includeInEachOverlapCheck, initialState.mutableData.moleculeIncluded);

            auto adjust = [&](const AngleSettings& settings,
                              const std::vector<GlycanShapePreference>& preferences,
                              const MutableData& mutableData,
                              const assembly::Selection& selection,
                              const std::vector<size_t>& glycanIds,
                              size_t times)
            {
                return adjustGlycans(
                    rng,
                    dihedralAngleDataTable,
                    overlapSettings,
                    settings,
                    times,
                    graph,
                    selection,
                    data,
                    preferences,
                    glycanIds,
                    mutableData);
            };
            uint cycle = 0;
            uint maxCycle = 0;
            GlycoproteinState bestState = initialState;
            bestState.totalOverlap = overlapVectorSum(
                totalOverlaps(overlapSettings, graph, data, fullSelection, bestState.mutableData.bounds));
            GlycoproteinState currentState = bestState;
            while (cycle < persistCycles && currentState.overlapSites.indices.size() > 0)
            {
                util::log(
                    __LINE__,
                    __FILE__,
                    util::INF,
                    "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles));
                AngleSettings settings = toAngleSettings(maxCycle);
                cycle++;
                maxCycle = std::max(cycle, maxCycle);
                std::vector<std::vector<size_t>> concerts = overlapConcerts(currentState.overlapSites);
                for (auto& concertGlycans : concerts)
                {
                    for (auto& glycanId : util::shuffleVector(rng, concertGlycans))
                    {
                        const std::vector<size_t>& linkageIds = data.glycans.linkages[glycanId];
                        currentState.preferences[glycanId] =
                            randomizeShape(rng, settings, data, currentState.mutableData, glycanId);
                        const GlycanShapePreference& glycanPreferences = currentState.preferences[glycanId];
                        for (size_t n = 0; n < linkageIds.size(); n++)
                        {
                            setLinkageShapeToPreference(
                                graph, data, currentState.mutableData, linkageIds[n], glycanPreferences[n]);
                        }
                    }

                    currentState.mutableData = adjust(
                        settings, currentState.preferences, currentState.mutableData, mainSelection, concertGlycans, 1);
                    adjustSidechains(rng, graph, data, currentState.mutableData, concertGlycans);
                    currentState.mutableData = adjust(
                        settings, currentState.preferences, currentState.mutableData, fullSelection, concertGlycans, 2);
                    currentState.totalOverlap = overlapVectorSum(
                        totalOverlaps(overlapSettings, graph, data, fullSelection, currentState.mutableData.bounds));
                    if (currentState.totalOverlap < bestState.totalOverlap)
                    {
                        cycle = 0;
                        bestState = currentState;
                    }
                    else
                    {
                        currentState = bestState;
                    }
                }
                currentState.overlapSites = determineOverlapState(
                    overlapSettings.rejectionThreshold,
                    overlapSettings,
                    graph,
                    data,
                    fullSelection,
                    currentState.mutableData.bounds);
            }
            return currentState;
        }

        GlycoproteinState resolveOverlapsWithWiggler(
            pcg32& rng,
            const DihedralAngleDataTable& dihedralAngleDataTable,
            PersistCycleAngleSettings toAngleSettings,
            SidechainAdjustment adjustSidechains,
            SidechainAdjustment restoreSidechains,
            GlycanShapeRandomizer& randomizeShape,
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const MutableData& initialState,
            size_t persistCycles,
            bool deleteSitesUntilResolved)
        {
            GlycoproteinState currentState;
            currentState.mutableData = initialState;
            AngleSettings initialAngleSettings = toAngleSettings(0);
            const std::vector<size_t> glycanIndices = util::indexVector(data.glycans.moleculeId);
            for (size_t glycanId : glycanIndices)
            {
                auto preference = randomizeShape(rng, initialAngleSettings, data, currentState.mutableData, glycanId);
                const std::vector<size_t>& linkageIds = data.glycans.linkages[glycanId];
                for (size_t k = 0; k < linkageIds.size(); k++)
                {
                    setLinkageShapeToPreference(graph, data, currentState.mutableData, linkageIds[k], preference[k]);
                }
                currentState.preferences.push_back(preference);
            }
            assembly::Selection fullSelection = assembly::selectByAtoms(graph, data.atoms.includeInEachOverlapCheck);
            currentState.overlapSites.indices = util::indexVector(data.glycans.moleculeId);
            currentState.mutableData = adjustGlycans(
                rng,
                dihedralAngleDataTable,
                overlapSettings,
                initialAngleSettings,
                2,
                graph,
                fullSelection,
                data,
                currentState.preferences,
                currentState.overlapSites.indices,
                currentState.mutableData);
            currentState.overlapSites = determineOverlapState(
                overlapSettings.rejectionThreshold,
                overlapSettings,
                graph,
                data,
                fullSelection,
                currentState.mutableData.bounds);
            for (bool done = false; !done;
                 done = currentState.overlapSites.indices.empty() || !deleteSitesUntilResolved)
            {
                currentState = randomDescent(
                    rng,
                    dihedralAngleDataTable,
                    toAngleSettings,
                    randomizeShape,
                    adjustSidechains,
                    persistCycles,
                    overlapSettings,
                    graph,
                    data,
                    currentState);
                if (deleteSitesUntilResolved && !currentState.overlapSites.indices.empty())
                {
                    size_t indexToRemove = util::randomIndex(rng, currentState.overlapSites.indices);
                    size_t glycan = currentState.overlapSites.indices[indexToRemove];
                    deleteMolecule(currentState.mutableData, data.glycans.moleculeId[glycan]);
                    fullSelection = assembly::selectByAtomsAndMolecules(
                        graph, data.atoms.includeInEachOverlapCheck, currentState.mutableData.moleculeIncluded);
                    size_t proteinResidue = data.glycans.attachmentResidue[glycan];
                    // restore atoms to initial shape
                    for (size_t n : residueAtoms(graph, proteinResidue))
                    {
                        currentState.mutableData.bounds.atoms[n] = data.atoms.initialState[n];
                    }
                    updateResidueBounds(graph, currentState.mutableData.bounds, proteinResidue);
                    updateResidueMoleculeBounds(graph, currentState.mutableData.bounds, proteinResidue);
                    currentState.overlapSites = determineOverlapState(
                        overlapSettings.rejectionThreshold,
                        overlapSettings,
                        graph,
                        data,
                        fullSelection,
                        currentState.mutableData.bounds);
                }
            }
            restoreSidechains(
                rng,
                graph,
                data,
                currentState.mutableData,
                includedGlycanIndices(data, currentState.mutableData.moleculeIncluded));
            util::log(__LINE__, __FILE__, util::INF, "Overlap: " + std::to_string(currentState.totalOverlap));
            return currentState;
        };
    } // namespace gpbuilder
} // namespace gmml
