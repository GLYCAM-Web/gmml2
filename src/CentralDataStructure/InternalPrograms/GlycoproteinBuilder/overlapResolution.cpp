#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapResolution.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinUtil.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/writerInterface.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapCount.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanWiggle.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/shapeRandomization.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/randomDescent.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/sidechains.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpBuilderStats.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileWriter.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileWriter.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/Assembly/assemblySelection.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/filesystem.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/structuredFiles.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/random.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/threads.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <variant>
#include <vector>
#include <string>
#include <sstream>
#include <ostream>
#include <stdexcept>

namespace glycoproteinBuilder
{
    void resolveOverlaps(const MolecularMetadata::SidechainRotamerData& sidechainRotamers,
                         const GlycoproteinAssembly& assembly, const GlycoproteinBuilderInputs& settings,
                         uint64_t rngSeed, const std::string& outputDir, const std::vector<std::string>& headerLines,
                         int numThreads)
    {
        const std::string structureDir = outputDir + "/samples";
        const std::string rejectDir    = structureDir + "/rejected";
        const std::string deletionDir  = structureDir + "/deletions";

        MetadataOrder initialMetadataOrder =
            [](pcg32&, const GlycamMetadata::DihedralAngleDataTable&, const std::vector<size_t>& metadataVector)
        {
            return codeUtils::indexVector(metadataVector);
        };
        MetadataOrder randomMetadataOrder = [](pcg32& rng, const GlycamMetadata::DihedralAngleDataTable& table,
                                               const std::vector<size_t>& metadataIndices)
        {
            return codeUtils::weightedRandomOrder(rng, codeUtils::indicesToValues(table.weights, metadataIndices));
        };
        bool freezeGlycositeResidueConformation = settings.useInitialGlycositeResidueConformation;
        auto standardDeviation = [](const AngleSettings& settings, const GlycamMetadata::DihedralAngleData& metadata)
        {
            std::function<std::pair<double, double>(const GlycamMetadata::AngleLimit&)> onLimit =
                [&](const GlycamMetadata::AngleLimit& dev)
            {
                double max_std   = settings.preferenceDeviation + settings.searchDeviation;
                double lower_std = dev.lowerDeviationLimit / max_std;
                double upper_std = dev.upperDeviationLimit / max_std;
                return std::pair<double, double> {lower_std, upper_std};
            };
            std::function<std::pair<double, double>(const GlycamMetadata::AngleStd&)> onStd =
                [&](const GlycamMetadata::AngleStd& dev)
            {
                return std::pair<double, double> {dev.lowerDeviationStd, dev.upperDeviationStd};
            };
            return cds::onAngleDeviation(onLimit, onStd, metadata.angle_deviation);
        };
        auto randomAngle = [&standardDeviation](pcg32& rng, const AngleSettings& settings,
                                                const GlycamMetadata::DihedralAngleData& metadata)
        {
            double stdCutoff = settings.preferenceDeviation;
            double num       = codeUtils::normalDistributionRandomDoubleWithCutoff(rng, -stdCutoff, stdCutoff);
            auto std         = standardDeviation(settings, metadata);
            return metadata.default_angle + num * (num < 0 ? std.first : std.second);
        };
        GlycanShapeRandomizer randomizeShape = [&randomAngle, &freezeGlycositeResidueConformation](
                                                   pcg32& rng, const AngleSettings& settings, const AssemblyData& data,
                                                   const MutableData& mutableData, size_t linkageId)
        {
            return randomLinkageShapePreference(rng, settings, data, mutableData.bounds, linkageId, randomAngle,
                                                freezeGlycositeResidueConformation);
        };

        auto getCoordinates = [](const std::vector<cds::Sphere>& bounds)
        {
            std::vector<cds::Coordinate> result;
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
                result.push_back(codeUtils::weightedRandomOrder(rng, data.residues.sidechainRotationWeights[residue]));
            }
            return result;
        };

        SidechainAdjustment adjustSidechains =
            [&](pcg32& rng, const AngleSettings& settings, WiggleGlycan wiggleGlycan, const assembly::Graph& graph,
                const AssemblyData& data, MutableData& mutableData,
                const std::vector<cds::GlycanShapePreference>& glycositePreferences, const std::vector<size_t>& glycans)
        {
            std::vector<size_t> sidechainResiduesWithGlycanOverlap;
            sidechainResiduesWithGlycanOverlap.reserve(graph.indices.residueCount);

            for (size_t residue : graph.residues.nodes.indices)
            {
                if (data.residues.types[residue] == cds::ResidueType::Protein &&
                    data.residues.sidechainDihedrals[residue].size() > 0 &&
                    sidechainHasGlycanOverlap(graph, data, mutableData, glycans, residue))
                {
                    sidechainResiduesWithGlycanOverlap.push_back(residue);
                }
            }
            std::vector<size_t> residues = codeUtils::shuffleVector(rng, sidechainResiduesWithGlycanOverlap);
            std::vector<std::vector<size_t>> preferences = sidechainRotationPreferences(rng, data, residues);
            for (size_t n = 0; n < residues.size(); n++)
            {
                setSidechainToLowestOverlapState(sidechainRotamers, graph, data, mutableData, preferences[n],
                                                 residues[n]);
            }
            assembly::Selection selection = assembly::selectByAtomsAndMolecules(
                graph, data.atoms.includeInEachOverlapCheck, mutableData.moleculeIncluded);
            for (size_t glycanId : codeUtils::shuffleVector(rng, glycans))
            {
                wiggleGlycan(graph, data, selection, settings, glycositePreferences[glycanId], mutableData, glycanId);
            }
        };

        SidechainAdjustment restoreSidechains =
            [&](pcg32& rng, const AngleSettings&, WiggleGlycan, const assembly::Graph& graph, const AssemblyData& data,
                MutableData& mutableData, const std::vector<cds::GlycanShapePreference>&,
                const std::vector<size_t>& glycans)
        {
            std::vector<size_t> residues = codeUtils::boolsToIndices(mutableData.residueSidechainMoved);
            for (size_t residue : residues)
            {
                restoreSidechainRotation(graph, data, mutableData, residue);
            }
            std::vector<std::vector<size_t>> preferences = sidechainRotationPreferences(rng, data, residues);
            for (size_t n = 0; n < residues.size(); n++)
            {
                if (sidechainHasGlycanOverlap(graph, data, mutableData, glycans, residues[n]))
                {
                    setSidechainToLowestOverlapState(sidechainRotamers, graph, data, mutableData, preferences[n],
                                                     residues[n]);
                }
            }
        };

        SidechainAdjustment noSidechainAdjustment =
            [&](pcg32&, const AngleSettings&, WiggleGlycan, const assembly::Graph&, const AssemblyData&, MutableData&,
                const std::vector<cds::GlycanShapePreference>&, const std::vector<size_t>&)
        {
            return;
        };

        WiggleGlycan wiggle = [&standardDeviation](const assembly::Graph& graph, const AssemblyData& data,
                                                   const assembly::Selection& selection, const AngleSettings& settings,
                                                   const cds::GlycanShapePreference& preference,
                                                   MutableData& mutableData, size_t glycanId)
        {
            auto searchAngles = [&standardDeviation, &settings](const GlycamMetadata::DihedralAngleData& metadata,
                                                                double preference, double deviation)
            {
                auto std = standardDeviation(settings, metadata);
                return cds::evenlySpacedAngles(preference, deviation * std.first, deviation * std.second,
                                               settings.searchIncrement);
            };
            return wiggleGlycan(graph, data, selection, {settings.searchDeviation, searchAngles}, preference,
                                mutableData, glycanId);
        };

        auto writeOffFile =
            [](const assembly::Graph& graph, const AssemblyData& data, const std::vector<cds::Coordinate>& coordinates,
               const std::vector<bool>& includedMolecules, const std::string& outputDir, const std::string& prefix)
        {
            std::vector<bool> residueIncluded(graph.indices.residueCount, false);
            for (size_t n = 0; n < graph.indices.residueCount; n++)
            {
                residueIncluded[n] = includedMolecules[graph.indices.residueMolecule[n]];
            }
            cds::OffFileData offData = toOffFileData(graph, data, coordinates);
            std::string fileName     = outputDir + "/" + prefix + ".off";
            codeUtils::writeToFile(fileName,
                                   [&](std::ostream& stream)
                                   {
                                       cds::WriteResiduesTogetherToOffFile(stream, graph, offData,
                                                                           codeUtils::boolsToIndices(residueIncluded),
                                                                           "GLYCOPROTEINBUILDER");
                                   });
        };

        auto residueTER = [](const AssemblyData& data, const std::vector<std::vector<size_t>>& residueIndices)
        {
            std::vector<std::vector<bool>> result;
            result.reserve(residueIndices.size());
            for (auto& indices : residueIndices)
            {
                result.push_back(cds::residueTER(codeUtils::indicesToValues(data.residues.types, indices)));
            }
            return result;
        };

        auto writePdbFile = [&](const assembly::Graph& graph, const AssemblyData& data,
                                const std::vector<cds::Coordinate>& coordinates, const std::vector<uint>& atomNumbers,
                                const std::vector<uint>& residueNumbers, const std::vector<bool>& includedMolecules,
                                const std::vector<std::array<size_t, 2>>& connectionIndices,
                                const std::string& outputDir, const std::string& prefix)
        {
            cds::PdbFileData pdbData = toPdbFileData(graph.indices, data, coordinates, atomNumbers, residueNumbers,
                                                     data.residues.chainIds, headerLines);
            std::vector<std::vector<size_t>> residueIndices =
                codeUtils::boolsToValues(graph.molecules.nodes.constituents, includedMolecules);
            std::vector<std::vector<bool>> TER = residueTER(data, residueIndices);
            codeUtils::createDirectories(outputDir);
            std::string filename = outputDir + "/" + prefix + ".pdb";
            codeUtils::writeToFile(filename,
                                   [&](std::ostream& stream)
                                   {
                                       cds::writeAssemblyToPdb(stream, graph, residueIndices, TER, connectionIndices,
                                                               pdbData);
                                       cds::theEnd(stream);
                                   });
        };

        const assembly::Graph& graph    = assembly.graph;
        const AssemblyData& data        = assembly.data;
        const MutableData& initialState = assembly.mutableData;
        SidechainAdjustment sidechainAdjustment =
            settings.moveOverlappingSidechains ? adjustSidechains : noSidechainAdjustment;
        SidechainAdjustment sidechainRestoration =
            settings.moveOverlappingSidechains ? restoreSidechains : noSidechainAdjustment;

        auto isNonProteinResidue = [&graph, &data](size_t n)
        {
            return data.molecules.types[graph.indices.residueMolecule[n]] != MoleculeType::protein;
        };

        std::vector<std::array<size_t, 2>> atomPairsConnectingNonProteinResidues;
        const graph::GraphEdges& residueEdges = graph.residues.edges;
        for (size_t n = 0; n < edgeCount(graph.residues); n++)
        {
            const std::array<size_t, 2>& adj = residueEdges.nodeAdjacencies[n];
            if (isNonProteinResidue(adj[0]) || isNonProteinResidue(adj[1]))
            {
                size_t atomEdgeId                   = residueEdges.indices[n];
                const std::array<size_t, 2> atomAdj = graph.atoms.edges.nodeAdjacencies[atomEdgeId];
                atomPairsConnectingNonProteinResidues.push_back(atomAdj);
            }
        }

        auto nonViableSites = [&settings](const assembly::Graph& graph, const assembly::Selection& selection,
                                          const AssemblyData& data, const std::vector<double>& overlaps)
        {
            std::vector<bool> result(graph.indices.moleculeCount, false);
            for (size_t molecule : includedGlycanMoleculeIds(data, selection.molecules))
            {
                for (size_t n : assembly::moleculeSelectedAtoms(graph, selection, molecule))
                {
                    if (cds::compareOverlaps(overlaps[n], settings.overlapRejectionThreshold) == 1)
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
            resolveOverlapsWithWiggler(rng, initialAngleSettings, mainAngleSettings, sidechainAdjustment,
                                       sidechainRestoration, randomizeShape, wiggle, graph, data, mutableData,
                                       settings.persistCycles, false);
            std::vector<cds::Coordinate> resolvedCoords = getCoordinates(mutableData.bounds.atoms);
            assembly::Selection selection               = assembly::selectAll(graph);
            std::vector<double> atomOverlaps            = totalOverlaps(graph, data, selection, mutableData.bounds);
            bool serialized                             = settings.MDprep;
            std::string prefix                          = "default";
            writePdbFile(graph, data, resolvedCoords, atomNumbers(serialized, data), residueNumbers(serialized, data),
                         mutableData.moleculeIncluded, atomPairsConnectingNonProteinResidues, outputDir, prefix);
            if (settings.writeOffFile)
            {
                writeOffFile(graph, data, resolvedCoords, mutableData.moleculeIncluded, outputDir, prefix);
            }
            std::vector<bool> nonViable = nonViableSites(graph, selection, data, atomOverlaps);
            stats                       = StructureStats {
                prefix,      false, false, mutableData.bounds, selection, nonViable, mutableData.residueSidechainMoved,
                atomOverlaps};
        };

        auto runIteration = [&](pcg32 rng, StructureStats& stats, const std::string& prefix)
        {
            MutableData mutableData = initialState;
            AngleSettings angleSettings {2.0, 0.5, 1.0, randomMetadataOrder};
            resolveOverlapsWithWiggler(rng, angleSettings, angleSettings, sidechainAdjustment, sidechainRestoration,
                                       randomizeShape, wiggle, graph, data, mutableData, settings.persistCycles,
                                       settings.deleteSitesUntilResolved);
            std::vector<cds::Coordinate> coordinates = getCoordinates(mutableData.bounds.atoms);
            const assembly::Selection fullSelection =
                assembly::selectByAtoms(graph, data.atoms.includeInEachOverlapCheck);

            bool hasDeleted               = codeUtils::contains(mutableData.moleculeIncluded, false);
            std::string directory         = outputDir;
            assembly::Selection selection = assembly::selectByAtomsAndMolecules(
                graph, data.atoms.includeInEachOverlapCheck, mutableData.moleculeIncluded);
            std::vector<double> atomOverlaps = totalOverlaps(graph, data, selection, mutableData.bounds);
            std::vector<bool> nonViable      = nonViableSites(graph, selection, data, atomOverlaps);
            bool reject                      = codeUtils::contains(nonViable, true);
            directory                        = hasDeleted ? deletionDir : (reject ? rejectDir : structureDir);
            bool serialized                  = settings.MDprep;
            writePdbFile(graph, data, coordinates, atomNumbers(serialized, data), residueNumbers(serialized, data),
                         mutableData.moleculeIncluded, atomPairsConnectingNonProteinResidues, directory, prefix);
            if (settings.writeOffFile)
            {
                writeOffFile(graph, data, coordinates, mutableData.moleculeIncluded, directory, prefix);
            }
            stats = StructureStats {prefix,
                                    reject,
                                    hasDeleted,
                                    mutableData.bounds,
                                    selection,
                                    nonViable,
                                    mutableData.residueSidechainMoved,
                                    atomOverlaps};
        };

        writePdbFile(graph, data, getCoordinates(data.atoms.initialState), data.atoms.numbers, data.residues.numbers,
                     initialState.moleculeIncluded, atomPairsConnectingNonProteinResidues, outputDir, "unresolved");

        // order of parallel for loop is undefined, so we generate all seeds up front for determinism
        size_t sampleCount     = settings.numberOfSamples;
        size_t totalStructures = sampleCount + 1;
        pcg32 seedingRng(rngSeed);
        std::vector<uint64_t> rngSeeds = codeUtils::randomIntegers<uint64_t>(totalStructures, seedingRng);
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
                std::string digits      = std::to_string(n);
                std::string extraZeroes = std::string(digitCount - digits.size(), '0');
                prefixes.push_back(extraZeroes + digits + "_glycoprotein");
            }
        }
        // will be initialized in loop
        std::vector<StructureStats> stats(totalStructures);

        codeUtils::setOpenMpNumberOfThreads(numThreads);
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

        Summary summary                                   = summarizeStats(graph, data, settings, rngSeed, stats);
        std::vector<codeUtils::TextVariant> textStructure = {
            codeUtils::TextVariant(codeUtils::TextHeader {2, "Glycoprotein Builder"}
             ),
            codeUtils::TextVariant(codeUtils::TextParagraph {headerLines}
             ),
            codeUtils::TextVariant(codeUtils::TextHeader {3, "Input"}
             ),
            codeUtils::TextVariant(
                codeUtils::TextParagraph {{"Filename: " + summary.filename, "Protein: " + summary.proteinFilename}}
             ),
            codeUtils::TextVariant(summary.parameterTable),
            codeUtils::TextVariant(codeUtils::TextHeader {3, "Structures"}
             ),
            codeUtils::TextVariant(summary.structuretable),
        };

        codeUtils::writeToFile(outputDir + "/summary.txt",
                               [&](std::ostream& stream)
                               {
                                   codeUtils::toTxt(stream, textStructure);
                               });
        codeUtils::writeToFile(outputDir + "/summary.html",
                               [&](std::ostream& stream)
                               {
                                   codeUtils::toHtml(stream, textStructure);
                               });
        codeUtils::writeToFile(outputDir + "/structures.csv",
                               [&](std::ostream& stream)
                               {
                                   codeUtils::toCsv(stream, ",", summary.structuretable);
                               });
    }
} // namespace glycoproteinBuilder
