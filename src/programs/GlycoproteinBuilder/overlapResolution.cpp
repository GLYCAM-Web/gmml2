#include "include/programs/GlycoproteinBuilder/overlapResolution.hpp"

#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyOverlap.hpp"
#include "include/assembly/assemblySelection.hpp"
#include "include/carbohydrate/dihedralAngleSearch.hpp"
#include "include/carbohydrate/dihedralShape.hpp"
#include "include/carbohydrate/pdbWriter.hpp"
#include "include/external/pcg/pcg_random.h"
#include "include/fileType/off/offFileData.hpp"
#include "include/fileType/off/offFileWriter.hpp"
#include "include/fileType/pdb/pdbFile.hpp"
#include "include/fileType/pdb/pdbFileData.hpp"
#include "include/fileType/pdb/pdbFileWriter.hpp"
#include "include/fileType/pdb/pdbResidue.hpp"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/glycoprotein/glycanShape.hpp"
#include "include/glycoprotein/glycanWiggle.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/glycoprotein/glycoproteinUtil.hpp"
#include "include/glycoprotein/overlapCount.hpp"
#include "include/glycoprotein/shapeRandomization.hpp"
#include "include/glycoprotein/sidechains.hpp"
#include "include/glycoprotein/writerInterface.hpp"
#include "include/graph/graphManipulation.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/programs/GlycoproteinBuilder/gpBuilderStats.hpp"
#include "include/programs/GlycoproteinBuilder/randomDescent.hpp"
#include "include/util/containers.hpp"
#include "include/util/files.hpp"
#include "include/util/filesystem.hpp"
#include "include/util/logging.hpp"
#include "include/util/random.hpp"
#include "include/util/strings.hpp"
#include "include/util/structuredFiles.hpp"
#include "include/util/threads.hpp"

#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        void resolveOverlaps(
            std::vector<StructureStats>& stats,
            const SidechainRotamerData& sidechainRotamers,
            const DihedralAngleDataTable& dihedralAngleDataTable,
            const GlycoproteinAssembly& assembly,
            const OverlapSettings& overlapSettings,
            const ResolutionSettings& resolutionSettings,
            const std::string& outputDir,
            const std::vector<std::string>& headerLines,
            int numThreads)
        {
            const std::string structureDir = outputDir + "/samples";
            const std::string rejectDir = structureDir + "/rejected";
            const std::string deletionDir = structureDir + "/deletions";

            MetadataOrder initialMetadataOrder =
                [](pcg32&, const DihedralAngleDataTable&, const std::vector<size_t>& metadataVector)
            { return util::indexVector(metadataVector); };
            MetadataOrder randomMetadataOrder =
                [](pcg32& rng, const DihedralAngleDataTable& table, const std::vector<size_t>& metadataIndices)
            { return util::weightedRandomOrder(rng, util::indicesToValues(table.weights, metadataIndices)); };
            bool freezeGlycositeResidueConformation = resolutionSettings.useInitialGlycositeResidueConformation;
            bool moveOverlappingSidechains = resolutionSettings.moveOverlappingSidechains;
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
            auto randomAngle =
                [&standardDeviation](pcg32& rng, const AngleSettings& settings, const DihedralAngleData& metadata)
            {
                double stdCutoff = settings.preferenceDeviation;
                double num = util::normalDistributionRandomDoubleWithCutoff(rng, -stdCutoff, stdCutoff);
                auto std = standardDeviation(settings, metadata);
                return metadata.default_angle + num * (num < 0 ? std.first : std.second);
            };
            GlycanShapeRandomizer randomizeShape =
                [&dihedralAngleDataTable, &randomAngle, &freezeGlycositeResidueConformation](
                    pcg32& rng,
                    const AngleSettings& settings,
                    const AssemblyData& data,
                    const MutableData& mutableData,
                    size_t linkageId)
            {
                return randomLinkageShapePreference(
                    rng,
                    dihedralAngleDataTable,
                    settings,
                    data,
                    mutableData.bounds,
                    linkageId,
                    randomAngle,
                    freezeGlycositeResidueConformation);
            };

            auto getCoordinates = [](const std::vector<Sphere>& bounds)
            {
                std::vector<Coordinate> result;
                result.reserve(bounds.size());
                for (size_t n = 0; n < bounds.size(); n++)
                {
                    result.push_back(bounds[n].center);
                }
                return result;
            };

            auto sidechainRotationPreferences =
                [](pcg32& rng, const AssemblyData& data, const std::vector<size_t>& residues)
            {
                std::vector<std::vector<size_t>> result;
                result.reserve(residues.size());
                for (size_t residue : residues)
                {
                    result.push_back(util::weightedRandomOrder(rng, data.residues.sidechainRotationWeights[residue]));
                }
                return result;
            };

            SidechainAdjustment adjustSidechains = [&](pcg32& rng,
                                                       const assembly::Graph& graph,
                                                       const AssemblyData& data,
                                                       MutableData& mutableData,
                                                       const std::vector<size_t>& glycans)
            {
                std::vector<size_t> sidechainResiduesWithGlycanOverlap;
                sidechainResiduesWithGlycanOverlap.reserve(residueCount(graph.source));

                for (size_t residue : indices(graph.residues.nodes))
                {
                    if (data.residues.types[residue] == ResidueType::Protein &&
                        data.residues.sidechainDihedrals[residue].size() > 0 &&
                        sidechainHasGlycanOverlap(overlapSettings, graph, data, mutableData, glycans, residue))
                    {
                        sidechainResiduesWithGlycanOverlap.push_back(residue);
                    }
                }
                std::vector<size_t> residues = util::shuffleVector(rng, sidechainResiduesWithGlycanOverlap);
                std::vector<std::vector<size_t>> preferences = sidechainRotationPreferences(rng, data, residues);
                for (size_t n = 0; n < residues.size(); n++)
                {
                    setSidechainToLowestOverlapState(
                        sidechainRotamers, overlapSettings, graph, data, mutableData, preferences[n], residues[n]);
                }
            };

            SidechainAdjustment restoreSidechains = [&](pcg32& rng,
                                                        const assembly::Graph& graph,
                                                        const AssemblyData& data,
                                                        MutableData& mutableData,
                                                        const std::vector<size_t>& glycans)
            {
                std::vector<size_t> residues = util::boolsToIndices(mutableData.residueSidechainMoved);
                for (size_t residue : residues)
                {
                    restoreSidechainRotation(graph, data, mutableData, residue);
                }
                std::vector<std::vector<size_t>> preferences = sidechainRotationPreferences(rng, data, residues);
                for (size_t n = 0; n < residues.size(); n++)
                {
                    if (sidechainHasGlycanOverlap(overlapSettings, graph, data, mutableData, glycans, residues[n]))
                    {
                        setSidechainToLowestOverlapState(
                            sidechainRotamers, overlapSettings, graph, data, mutableData, preferences[n], residues[n]);
                    }
                }
            };

            SidechainAdjustment noSidechainAdjustment =
                [&](pcg32&, const assembly::Graph&, const AssemblyData&, MutableData&, const std::vector<size_t>&)
            { return; };

            auto writeOffFile = [](const assembly::Graph& graph,
                                   const AssemblyData& data,
                                   const std::vector<Coordinate>& coordinates,
                                   const std::vector<bool>& includedMolecules,
                                   const std::string& outputDir,
                                   const std::string& prefix)
            {
                std::vector<bool> residueIncluded(residueCount(graph.source), false);
                for (size_t n = 0; n < residueCount(graph.source); n++)
                {
                    residueIncluded[n] = includedMolecules[residueMolecule(graph.source.indices, n)];
                }
                off::OffFileData offData = toOffFileData(graph, data, coordinates);
                std::string fileName = outputDir + "/" + prefix + ".off";
                util::writeToFile(
                    fileName,
                    [&](std::ostream& stream) {
                        off::writeResiduesTogether(
                            stream, offData, graph, util::boolsToIndices(residueIncluded), "GLYCOPROTEINBUILDER");
                    });
            };

            auto residueTER = [](const AssemblyData& data, const std::vector<std::vector<size_t>>& residueIndices)
            {
                std::vector<std::vector<bool>> result;
                result.reserve(residueIndices.size());
                for (auto& indices : residueIndices)
                {
                    result.push_back(pdb::residueTER(util::indicesToValues(data.residues.types, indices)));
                }
                return result;
            };

            auto writePdbFile = [&](const assembly::Graph& graph,
                                    const AssemblyData& data,
                                    const std::vector<Coordinate>& coordinates,
                                    const std::vector<uint>& atomNumbers,
                                    const std::vector<uint>& residueNumbers,
                                    const std::vector<bool>& includedMolecules,
                                    const std::vector<std::array<size_t, 2>>& connectionIndices,
                                    const std::string& outputDir,
                                    const std::string& prefix)
            {
                pdb::PdbFileData pdbData = toPdbFileData(
                    graph.source.indices,
                    data,
                    coordinates,
                    atomNumbers,
                    residueNumbers,
                    data.residues.chainIds,
                    headerLines);
                std::vector<std::vector<size_t>> residueIndices =
                    util::boolsToValues(graph.molecules.nodes.constituents, includedMolecules);
                std::vector<std::vector<bool>> TER = residueTER(data, residueIndices);
                util::createDirectories(outputDir);
                std::string filename = outputDir + "/" + prefix + ".pdb";
                util::writeToFile(
                    filename,
                    [&](std::ostream& stream)
                    {
                        pdb::writeAssemblyToPdb(stream, pdbData, graph, residueIndices, TER, connectionIndices);
                        pdb::theEnd(stream);
                    });
            };

            const assembly::Graph& graph = assembly.graph;
            const AssemblyData& data = assembly.data;
            const MutableData& initialState = assembly.mutableData;
            SidechainAdjustment sidechainAdjustment =
                moveOverlappingSidechains ? adjustSidechains : noSidechainAdjustment;
            SidechainAdjustment sidechainRestoration =
                moveOverlappingSidechains ? restoreSidechains : noSidechainAdjustment;

            auto isNonProteinResidue = [&graph, &data](size_t n)
            { return data.molecules.types[residueMolecule(graph.source.indices, n)] != MoleculeType::protein; };

            std::vector<std::array<size_t, 2>> atomPairsConnectingNonProteinResidues;
            for (size_t n = 0; n < edgeCount(graph.residues); n++)
            {
                const std::array<size_t, 2>& adj = graph.residues.edges.nodeAdjacencies[n];
                if (isNonProteinResidue(adj[0]) || isNonProteinResidue(adj[1]))
                {
                    size_t atomEdgeId = sourceIndex(graph.residues.edges, n);
                    const std::array<size_t, 2> atomAdj = graph.atoms.edges.nodeAdjacencies[atomEdgeId];
                    atomPairsConnectingNonProteinResidues.push_back(atomAdj);
                }
            }

            auto nonViableSites = [&overlapSettings](
                                      const assembly::Graph& graph,
                                      const assembly::Selection& selection,
                                      const AssemblyData& data,
                                      const std::vector<double>& overlaps)
            {
                std::vector<bool> result(moleculeCount(graph.source), false);
                for (size_t molecule : includedGlycanMoleculeIds(data, selection.molecules))
                {
                    for (size_t n : assembly::moleculeSelectedAtoms(graph, selection, molecule))
                    {
                        if (compareOverlaps(overlaps[n], overlapSettings.rejectionThreshold) == 1)
                        {
                            result[molecule] = true;
                        }
                    }
                }
                return result;
            };

            auto runInitial = [&](pcg32 rng, StructureStats& stats)
            {
                AngleSettings initialAngleSettings {0.0, 2.0, 1.0, initialMetadataOrder};
                AngleSettings mainAngleSettings {0.5, 2.0, 1.0, initialMetadataOrder};
                GlycoproteinState state = resolveOverlapsWithWiggler(
                    rng,
                    dihedralAngleDataTable,
                    initialAngleSettings,
                    mainAngleSettings,
                    sidechainAdjustment,
                    sidechainRestoration,
                    randomizeShape,
                    overlapSettings,
                    graph,
                    data,
                    initialState,
                    resolutionSettings.persistCycles,
                    false);
                std::vector<Coordinate> resolvedCoords = getCoordinates(state.mutableData.bounds.atoms);
                assembly::Selection selection = assembly::selectAll(graph);
                std::vector<double> atomOverlaps =
                    totalOverlaps(overlapSettings, graph, data, selection, state.mutableData.bounds);
                bool serialized = resolutionSettings.prepareForMD;
                std::string prefix = "default";
                writePdbFile(
                    graph,
                    data,
                    resolvedCoords,
                    atomNumbers(serialized, data),
                    residueNumbers(serialized, data),
                    state.mutableData.moleculeIncluded,
                    atomPairsConnectingNonProteinResidues,
                    outputDir,
                    prefix);
                if (resolutionSettings.writeOffFile)
                {
                    writeOffFile(graph, data, resolvedCoords, state.mutableData.moleculeIncluded, outputDir, prefix);
                }
                std::vector<bool> nonViable = nonViableSites(graph, selection, data, atomOverlaps);
                stats = StructureStats {
                    prefix,
                    false,
                    false,
                    state.mutableData.bounds,
                    selection,
                    nonViable,
                    state.mutableData.residueSidechainMoved,
                    atomOverlaps};
            };

            auto runIteration = [&](pcg32 rng, StructureStats& stats, const std::string& prefix)
            {
                AngleSettings angleSettings {2.0, 0.5, 1.0, randomMetadataOrder};
                GlycoproteinState state = resolveOverlapsWithWiggler(
                    rng,
                    dihedralAngleDataTable,
                    angleSettings,
                    angleSettings,
                    sidechainAdjustment,
                    sidechainRestoration,
                    randomizeShape,
                    overlapSettings,
                    graph,
                    data,
                    initialState,
                    resolutionSettings.persistCycles,
                    resolutionSettings.deleteSitesUntilResolved);
                std::vector<Coordinate> coordinates = getCoordinates(state.mutableData.bounds.atoms);
                const assembly::Selection fullSelection =
                    assembly::selectByAtoms(graph, data.atoms.includeInEachOverlapCheck);

                bool hasDeleted = util::contains(state.mutableData.moleculeIncluded, false);
                std::string directory = outputDir;
                assembly::Selection selection = assembly::selectByAtomsAndMolecules(
                    graph, data.atoms.includeInEachOverlapCheck, state.mutableData.moleculeIncluded);
                std::vector<double> atomOverlaps =
                    totalOverlaps(overlapSettings, graph, data, selection, state.mutableData.bounds);
                std::vector<bool> nonViable = nonViableSites(graph, selection, data, atomOverlaps);
                bool reject = util::contains(nonViable, true);
                directory = hasDeleted ? deletionDir : (reject ? rejectDir : structureDir);
                bool serialized = resolutionSettings.prepareForMD;
                writePdbFile(
                    graph,
                    data,
                    coordinates,
                    atomNumbers(serialized, data),
                    residueNumbers(serialized, data),
                    state.mutableData.moleculeIncluded,
                    atomPairsConnectingNonProteinResidues,
                    directory,
                    prefix);
                if (resolutionSettings.writeOffFile)
                {
                    writeOffFile(graph, data, coordinates, state.mutableData.moleculeIncluded, directory, prefix);
                }
                stats = StructureStats {
                    prefix,
                    reject,
                    hasDeleted,
                    state.mutableData.bounds,
                    selection,
                    nonViable,
                    state.mutableData.residueSidechainMoved,
                    atomOverlaps};
            };

            writePdbFile(
                graph,
                data,
                getCoordinates(data.atoms.initialState),
                data.atoms.numbers,
                data.residues.numbers,
                initialState.moleculeIncluded,
                atomPairsConnectingNonProteinResidues,
                outputDir,
                "unresolved");

            // order of parallel for loop is undefined, so we generate all seeds up front for determinism
            size_t sampleCount = resolutionSettings.numberOfSamples;
            size_t totalStructures = sampleCount + 1;
            pcg32 seedingRng(resolutionSettings.rngSeed);
            std::vector<uint64_t> rngSeeds = util::randomIntegers<uint64_t>(totalStructures, seedingRng);
            std::vector<std::string> prefixes;
            prefixes.reserve(sampleCount);
            if (sampleCount == 1)
            {
                prefixes.push_back("glycoprotein");
            }
            else if (sampleCount > 1)
            {
                size_t digitCount = std::to_string(sampleCount - 1).size();
                for (size_t n = 0; n < sampleCount; n++)
                {
                    std::string digits = std::to_string(n);
                    std::string extraZeroes = std::string(digitCount - digits.size(), '0');
                    prefixes.push_back(extraZeroes + digits + "_glycoprotein");
                }
            }
            // will be initialized in loop
            stats.resize(totalStructures);

            util::setOpenMpNumberOfThreads(numThreads);
//        clang format doesn't align pragmas
/*        */#pragma omp parallel for
            for (size_t count = 0; count < totalStructures; count++)
            {
                pcg32 rng(rngSeeds[count]);
                if (count == 0)
                {
                    runInitial(rng, stats[0]);
                }
                else
                {
                    runIteration(rng, stats[count], prefixes[count - 1]);
                }
            }
        }
    } // namespace gpbuilder
} // namespace gmml
