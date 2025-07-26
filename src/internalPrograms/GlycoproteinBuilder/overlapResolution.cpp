#include "include/internalPrograms/GlycoproteinBuilder/overlapResolution.hpp"

#include "include/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "include/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/CentralDataStructure/pdbWriter.hpp"
#include "include/External_Libraries/PCG/pcg_random.h"
#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblySelection.hpp"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/graph/graphManipulation.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycanWiggle.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinUtil.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/gpBuilderStats.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/overlapCount.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/randomDescent.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/shapeRandomization.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/sidechains.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/writerInterface.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/off/offFileData.hpp"
#include "include/off/offFileWriter.hpp"
#include "include/pdb/pdbFile.hpp"
#include "include/pdb/pdbFileData.hpp"
#include "include/pdb/pdbFileWriter.hpp"
#include "include/pdb/pdbResidue.hpp"
#include "include/structure/atomOverlaps.hpp"
#include "include/util/casting.hpp"
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
            const SidechainRotamerData& sidechainRotamers,
            const GlycoproteinAssembly& assembly,
            const GlycoproteinBuilderInputs& settings,
            uint64_t rngSeed,
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
            bool freezeGlycositeResidueConformation = settings.useInitialGlycositeResidueConformation;
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
            GlycanShapeRandomizer randomizeShape = [&randomAngle, &freezeGlycositeResidueConformation](
                                                       pcg32& rng,
                                                       const AngleSettings& settings,
                                                       const AssemblyData& data,
                                                       const MutableData& mutableData,
                                                       size_t linkageId)
            {
                return randomLinkageShapePreference(
                    rng,
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
                                                       const AngleSettings& settings,
                                                       WiggleGlycan wiggleGlycan,
                                                       const assembly::Graph& graph,
                                                       const AssemblyData& data,
                                                       MutableData& mutableData,
                                                       const std::vector<GlycanShapePreference>& glycositePreferences,
                                                       const std::vector<size_t>& glycans)
            {
                std::vector<size_t> sidechainResiduesWithGlycanOverlap;
                sidechainResiduesWithGlycanOverlap.reserve(graph.indices.residueCount);

                for (size_t residue : graph.residues.nodes.indices)
                {
                    if (data.residues.types[residue] == ResidueType::Protein &&
                        data.residues.sidechainDihedrals[residue].size() > 0 &&
                        sidechainHasGlycanOverlap(graph, data, mutableData, glycans, residue))
                    {
                        sidechainResiduesWithGlycanOverlap.push_back(residue);
                    }
                }
                std::vector<size_t> residues = util::shuffleVector(rng, sidechainResiduesWithGlycanOverlap);
                std::vector<std::vector<size_t>> preferences = sidechainRotationPreferences(rng, data, residues);
                for (size_t n = 0; n < residues.size(); n++)
                {
                    setSidechainToLowestOverlapState(
                        sidechainRotamers, graph, data, mutableData, preferences[n], residues[n]);
                }
                assembly::Selection selection = assembly::selectByAtomsAndMolecules(
                    graph, data.atoms.includeInEachOverlapCheck, mutableData.moleculeIncluded);
                for (size_t glycanId : util::shuffleVector(rng, glycans))
                {
                    wiggleGlycan(
                        graph, data, selection, settings, glycositePreferences[glycanId], mutableData, glycanId);
                }
            };

            SidechainAdjustment restoreSidechains = [&](pcg32& rng,
                                                        const AngleSettings&,
                                                        WiggleGlycan,
                                                        const assembly::Graph& graph,
                                                        const AssemblyData& data,
                                                        MutableData& mutableData,
                                                        const std::vector<GlycanShapePreference>&,
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
                    if (sidechainHasGlycanOverlap(graph, data, mutableData, glycans, residues[n]))
                    {
                        setSidechainToLowestOverlapState(
                            sidechainRotamers, graph, data, mutableData, preferences[n], residues[n]);
                    }
                }
            };

            SidechainAdjustment noSidechainAdjustment = [&](pcg32&,
                                                            const AngleSettings&,
                                                            WiggleGlycan,
                                                            const assembly::Graph&,
                                                            const AssemblyData&,
                                                            MutableData&,
                                                            const std::vector<GlycanShapePreference>&,
                                                            const std::vector<size_t>&) { return; };

            WiggleGlycan wiggle = [&standardDeviation](
                                      const assembly::Graph& graph,
                                      const AssemblyData& data,
                                      const assembly::Selection& selection,
                                      const AngleSettings& settings,
                                      const GlycanShapePreference& preference,
                                      MutableData& mutableData,
                                      size_t glycanId)
            {
                auto searchAngles = [&standardDeviation,
                                     &settings](const DihedralAngleData& metadata, double preference, double deviation)
                {
                    auto std = standardDeviation(settings, metadata);
                    return evenlySpacedAngles(
                        preference, deviation * std.first, deviation * std.second, settings.searchIncrement);
                };
                return wiggleGlycan(
                    graph,
                    data,
                    selection,
                    {settings.searchDeviation, searchAngles},
                    preference,
                    mutableData,
                    glycanId);
            };

            auto writeOffFile = [](const assembly::Graph& graph,
                                   const AssemblyData& data,
                                   const std::vector<Coordinate>& coordinates,
                                   const std::vector<bool>& includedMolecules,
                                   const std::string& outputDir,
                                   const std::string& prefix)
            {
                std::vector<bool> residueIncluded(graph.indices.residueCount, false);
                for (size_t n = 0; n < graph.indices.residueCount; n++)
                {
                    residueIncluded[n] = includedMolecules[graph.indices.residueMolecule[n]];
                }
                off::OffFileData offData = toOffFileData(graph, data, coordinates);
                std::string fileName = outputDir + "/" + prefix + ".off";
                util::writeToFile(
                    fileName,
                    [&](std::ostream& stream) {
                        off::writeResiduesTogether(
                            stream, graph, offData, util::boolsToIndices(residueIncluded), "GLYCOPROTEINBUILDER");
                    });
            };

            auto residueTER = [](const AssemblyData& data, const std::vector<std::vector<size_t>>& residueIndices)
            {
                std::vector<std::vector<bool>> result;
                result.reserve(residueIndices.size());
                for (auto& indices : residueIndices)
                {
                    result.push_back(gmml::residueTER(util::indicesToValues(data.residues.types, indices)));
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
                    graph.indices, data, coordinates, atomNumbers, residueNumbers, data.residues.chainIds, headerLines);
                std::vector<std::vector<size_t>> residueIndices =
                    util::boolsToValues(graph.molecules.nodes.constituents, includedMolecules);
                std::vector<std::vector<bool>> TER = residueTER(data, residueIndices);
                util::createDirectories(outputDir);
                std::string filename = outputDir + "/" + prefix + ".pdb";
                util::writeToFile(
                    filename,
                    [&](std::ostream& stream)
                    {
                        pdb::writeAssemblyToPdb(stream, graph, residueIndices, TER, connectionIndices, pdbData);
                        pdb::theEnd(stream);
                    });
            };

            const assembly::Graph& graph = assembly.graph;
            const AssemblyData& data = assembly.data;
            const MutableData& initialState = assembly.mutableData;
            SidechainAdjustment sidechainAdjustment =
                settings.moveOverlappingSidechains ? adjustSidechains : noSidechainAdjustment;
            SidechainAdjustment sidechainRestoration =
                settings.moveOverlappingSidechains ? restoreSidechains : noSidechainAdjustment;

            auto isNonProteinResidue = [&graph, &data](size_t n)
            { return data.molecules.types[graph.indices.residueMolecule[n]] != MoleculeType::protein; };

            std::vector<std::array<size_t, 2>> atomPairsConnectingNonProteinResidues;
            const graph::GraphEdges& residueEdges = graph.residues.edges;
            for (size_t n = 0; n < edgeCount(graph.residues); n++)
            {
                const std::array<size_t, 2>& adj = residueEdges.nodeAdjacencies[n];
                if (isNonProteinResidue(adj[0]) || isNonProteinResidue(adj[1]))
                {
                    size_t atomEdgeId = residueEdges.indices[n];
                    const std::array<size_t, 2> atomAdj = graph.atoms.edges.nodeAdjacencies[atomEdgeId];
                    atomPairsConnectingNonProteinResidues.push_back(atomAdj);
                }
            }

            auto nonViableSites = [&settings](
                                      const assembly::Graph& graph,
                                      const assembly::Selection& selection,
                                      const AssemblyData& data,
                                      const std::vector<double>& overlaps)
            {
                std::vector<bool> result(graph.indices.moleculeCount, false);
                for (size_t molecule : includedGlycanMoleculeIds(data, selection.molecules))
                {
                    for (size_t n : assembly::moleculeSelectedAtoms(graph, selection, molecule))
                    {
                        if (compareOverlaps(overlaps[n], settings.overlapRejectionThreshold) == 1)
                        {
                            result[molecule] = true;
                        }
                    }
                }
                return result;
            };

            auto runInitial = [&](pcg32 rng, StructureStats& stats)
            {
                MutableData mutableData = initialState;
                AngleSettings initialAngleSettings {0.0, 2.0, 1.0, initialMetadataOrder};
                AngleSettings mainAngleSettings {0.5, 2.0, 1.0, initialMetadataOrder};
                resolveOverlapsWithWiggler(
                    rng,
                    initialAngleSettings,
                    mainAngleSettings,
                    sidechainAdjustment,
                    sidechainRestoration,
                    randomizeShape,
                    wiggle,
                    graph,
                    data,
                    mutableData,
                    settings.persistCycles,
                    false);
                std::vector<Coordinate> resolvedCoords = getCoordinates(mutableData.bounds.atoms);
                assembly::Selection selection = assembly::selectAll(graph);
                std::vector<double> atomOverlaps = totalOverlaps(graph, data, selection, mutableData.bounds);
                bool serialized = settings.MDprep;
                std::string prefix = "default";
                writePdbFile(
                    graph,
                    data,
                    resolvedCoords,
                    atomNumbers(serialized, data),
                    residueNumbers(serialized, data),
                    mutableData.moleculeIncluded,
                    atomPairsConnectingNonProteinResidues,
                    outputDir,
                    prefix);
                if (settings.writeOffFile)
                {
                    writeOffFile(graph, data, resolvedCoords, mutableData.moleculeIncluded, outputDir, prefix);
                }
                std::vector<bool> nonViable = nonViableSites(graph, selection, data, atomOverlaps);
                stats = StructureStats {
                    prefix,
                    false,
                    false,
                    mutableData.bounds,
                    selection,
                    nonViable,
                    mutableData.residueSidechainMoved,
                    atomOverlaps};
            };

            auto runIteration = [&](pcg32 rng, StructureStats& stats, const std::string& prefix)
            {
                MutableData mutableData = initialState;
                AngleSettings angleSettings {2.0, 0.5, 1.0, randomMetadataOrder};
                resolveOverlapsWithWiggler(
                    rng,
                    angleSettings,
                    angleSettings,
                    sidechainAdjustment,
                    sidechainRestoration,
                    randomizeShape,
                    wiggle,
                    graph,
                    data,
                    mutableData,
                    settings.persistCycles,
                    settings.deleteSitesUntilResolved);
                std::vector<Coordinate> coordinates = getCoordinates(mutableData.bounds.atoms);
                const assembly::Selection fullSelection =
                    assembly::selectByAtoms(graph, data.atoms.includeInEachOverlapCheck);

                bool hasDeleted = util::contains(mutableData.moleculeIncluded, false);
                std::string directory = outputDir;
                assembly::Selection selection = assembly::selectByAtomsAndMolecules(
                    graph, data.atoms.includeInEachOverlapCheck, mutableData.moleculeIncluded);
                std::vector<double> atomOverlaps = totalOverlaps(graph, data, selection, mutableData.bounds);
                std::vector<bool> nonViable = nonViableSites(graph, selection, data, atomOverlaps);
                bool reject = util::contains(nonViable, true);
                directory = hasDeleted ? deletionDir : (reject ? rejectDir : structureDir);
                bool serialized = settings.MDprep;
                writePdbFile(
                    graph,
                    data,
                    coordinates,
                    atomNumbers(serialized, data),
                    residueNumbers(serialized, data),
                    mutableData.moleculeIncluded,
                    atomPairsConnectingNonProteinResidues,
                    directory,
                    prefix);
                if (settings.writeOffFile)
                {
                    writeOffFile(graph, data, coordinates, mutableData.moleculeIncluded, directory, prefix);
                }
                stats = StructureStats {
                    prefix,
                    reject,
                    hasDeleted,
                    mutableData.bounds,
                    selection,
                    nonViable,
                    mutableData.residueSidechainMoved,
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
            size_t sampleCount = settings.numberOfSamples;
            size_t totalStructures = sampleCount + 1;
            pcg32 seedingRng(rngSeed);
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
            std::vector<StructureStats> stats(totalStructures);

            util::setOpenMpNumberOfThreads(numThreads);
// clang format doesn't align pragmas
/*    */#pragma omp parallel for
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

            Summary summary = summarizeStats(graph, data, settings, rngSeed, stats);
            std::vector<util::TextVariant> textStructure = {
                util::TextVariant(util::TextHeader {2, "Glycoprotein Builder"}
                 ),
                util::TextVariant(util::TextParagraph {headerLines}
                 ),
                util::TextVariant(util::TextHeader {3, "Input"}
                 ),
                util::TextVariant(
                    util::TextParagraph {{"Filename: " + summary.filename, "Protein: " + summary.proteinFilename}}
                 ),
                util::TextVariant(summary.parameterTable),
                util::TextVariant(util::TextHeader {3, "Structures"}
                 ),
                util::TextVariant(summary.structuretable),
            };

            util::writeToFile(
                outputDir + "/summary.txt", [&](std::ostream& stream) { util::toTxt(stream, textStructure); });
            util::writeToFile(
                outputDir + "/summary.html", [&](std::ostream& stream) { util::toHtml(stream, textStructure); });
            util::writeToFile(
                outputDir + "/structures.csv",
                [&](std::ostream& stream) { util::toCsv(stream, ",", summary.structuretable); });
        }
    } // namespace gpbuilder
} // namespace gmml
