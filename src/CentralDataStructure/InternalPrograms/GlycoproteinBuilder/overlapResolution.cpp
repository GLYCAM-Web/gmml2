#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapResolution.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/writerInterface.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapCount.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/shapeRandomization.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/randomDescent.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/sidechains.hpp"
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
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/files.hpp"
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
#include <fstream>
#include <stdexcept>

namespace glycoproteinBuilder
{
    void resolveOverlaps(const MolecularMetadata::SidechainRotamerData& sidechainRotamers,
                         const OverlapMultiplier& overlapMultiplier, const GlycoproteinAssembly& assembly,
                         const GlycoproteinBuilderInputs& settings, const std::string& outputDir,
                         const std::vector<std::string>& headerLines, int numThreads)
    {
        const std::string rejectDir   = outputDir + "/rejected";
        const std::string deletionDir = outputDir + "/deletions";
        uint64_t seed                 = settings.isDeterministic ? settings.seed : codeUtils::generateRandomSeed();
        pcg32 seedingRng(seed);

        auto randomMetadata = [](pcg32& rng, GlycamMetadata::DihedralAngleDataVector metadataVector)
        {
            auto weights = GlycamMetadata::dihedralAngleDataWeights(metadataVector);
            return codeUtils::weightedRandomOrder(rng, weights);
        };
        bool freezeGlycositeResidueConformation = settings.freezeGlycositeResidueConformation;
        double preferenceDeviation              = 2.0;
        double searchDeviation                  = 0.5;
        double searchIncrement                  = 1.0;
        auto standardDeviation =
            [&preferenceDeviation, &searchDeviation](const GlycamMetadata::DihedralAngleData& metadata)
        {
            std::function<std::pair<double, double>(const GlycamMetadata::AngleLimit&)> onLimit =
                [&](const GlycamMetadata::AngleLimit& dev)
            {
                double max_std   = preferenceDeviation + searchDeviation;
                double lower_std = dev.lowerDeviationLimit / max_std;
                double upper_std = dev.upperDeviationLimit / max_std;
                return std::pair<double, double> {lower_std, upper_std};
            };
            std::function<std::pair<double, double>(const GlycamMetadata::AngleStd&)> onStd =
                [&](const GlycamMetadata::AngleStd& dev)
            {
                return std::pair<double, double> {dev.lowerDeviationStd, dev.upperDeviationStd};
            };
            return GlycamMetadata::onAngleDeviation(onLimit, onStd, metadata.angle_deviation);
        };
        auto searchAngles = [&standardDeviation, &searchIncrement](const GlycamMetadata::DihedralAngleData& metadata,
                                                                   double preference, double deviation)
        {
            auto std = standardDeviation(metadata);
            return cds::evenlySpacedAngles(preference, deviation * std.first, deviation * std.second, searchIncrement);
        };
        cds::AngleSearchSettings searchSettings = cds::AngleSearchSettings {searchDeviation, searchAngles};
        auto randomAngle =
            [&standardDeviation, &preferenceDeviation](pcg32& rng, GlycamMetadata::DihedralAngleData metadata)
        {
            double stdCutoff = preferenceDeviation;
            double num       = codeUtils::normalDistributionRandomDoubleWithCutoff(rng, -stdCutoff, stdCutoff);
            auto std         = standardDeviation(metadata);
            return metadata.default_angle + num * (num < 0 ? std.first : std.second);
        };
        auto randomizeShape = [&randomMetadata, &randomAngle,
                               &freezeGlycositeResidueConformation](pcg32& rng, const AssemblyData& data,
                                                                    const MutableData& mutableData, size_t linkageId)
        {
            return randomLinkageShapePreference(rng, data, mutableData.bounds, linkageId, randomMetadata, randomAngle,
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

        auto deleteMolecule = [](MutableData& mutableData, size_t moleculeId)
        {
            mutableData.moleculeIncluded[moleculeId] = false;
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
            [&](pcg32& rng, const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                const std::vector<cds::GlycanShapePreference>& glycositePreferences, const std::vector<size_t>& glycans)
        {
            std::vector<size_t> sidechainResiduesWithGlycanOverlap;
            sidechainResiduesWithGlycanOverlap.reserve(graph.residueCount);

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
            for (size_t glycanId : codeUtils::shuffleVector(rng, glycans))
            {
                wiggleGlycan(graph, data, mutableData, data.atoms.includeInEachOverlapCheck, glycanId, searchSettings,
                             overlapMultiplier, glycositePreferences[glycanId]);
            }
        };

        SidechainAdjustment restoreSidechains =
            [&](pcg32& rng, const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                const std::vector<cds::GlycanShapePreference>&, const std::vector<size_t>& glycans)
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

        SidechainAdjustment noSidechainAdjustment = [&](pcg32&, const assembly::Graph&, const AssemblyData&,
                                                        MutableData&, const std::vector<cds::GlycanShapePreference>&,
                                                        const std::vector<size_t>&)
        {
            return;
        };

        auto resolveOverlapsWithWiggler = [&](pcg32& rng, SidechainAdjustment adjustSidechains,
                                              SidechainAdjustment restoreSidechains, const assembly::Graph& graph,
                                              const AssemblyData& data, MutableData& mutableData,
                                              bool deleteSitesUntilResolved)
        {
            std::vector<cds::GlycanShapePreference> glycositePreferences;
            const std::vector<size_t> glycanIndices = codeUtils::indexVector(data.glycans.moleculeId);
            for (size_t glycanId : glycanIndices)
            {
                auto preference                       = randomizeShape(rng, data, mutableData, glycanId);
                const std::vector<size_t>& linkageIds = data.glycans.linkages[glycanId];
                for (size_t k = 0; k < linkageIds.size(); k++)
                {
                    setLinkageShapeToPreference(graph, data, mutableData, linkageIds[k], preference[k]);
                }
                updateGlycanBounds(graph, data, mutableData.bounds, glycanId);
                glycositePreferences.push_back(preference);
            }
            for (size_t glycanId : codeUtils::shuffleVector(rng, glycanIndices))
            {
                wiggleGlycan(graph, data, mutableData, data.atoms.includeInMainOverlapCheck, glycanId, searchSettings,
                             overlapMultiplier, glycositePreferences[glycanId]);
            }
            adjustSidechains(rng, graph, data, mutableData, glycositePreferences, glycanIndices);
            GlycoproteinState currentState;
            std::vector<size_t> overlapSites = determineSitesWithOverlap(glycanIndices, graph, data, mutableData,
                                                                         data.atoms.includeInEachOverlapCheck);
            for (bool done = false; !done; done = overlapSites.empty() || !deleteSitesUntilResolved)
            {
                cds::Overlap initialOverlap =
                    cds::overlapVectorSum(totalOverlaps(graph, data, mutableData, data.defaultResidueWeight,
                                                        data.atoms.includeInEachOverlapCheck, overlapMultiplier));

                GlycoproteinState initialState = {initialOverlap, overlapSites, glycositePreferences};
                currentState =
                    randomDescent(rng, randomizeShape, adjustSidechains, searchSettings, settings.persistCycles,
                                  overlapMultiplier, graph, data, mutableData, initialState);
                overlapSites = currentState.overlapSites;
                if (deleteSitesUntilResolved && !overlapSites.empty())
                {
                    size_t indexToRemove = codeUtils::randomIndex(rng, overlapSites);
                    size_t glycan        = overlapSites[indexToRemove];
                    deleteMolecule(mutableData, data.glycans.moleculeId[glycan]);
                    size_t proteinResidue = data.glycans.attachmentResidue[glycan];
                    // restore atoms to initial shape
                    for (size_t n : residueAtoms(graph, proteinResidue))
                    {
                        mutableData.bounds.atoms[n] = data.atoms.initialState[n];
                    }
                    updateResidueBounds(graph, mutableData.bounds, proteinResidue);
                    updateResidueMoleculeBounds(graph, mutableData.bounds, proteinResidue);
                    overlapSites = determineSitesWithOverlap(includedGlycanIndices(data, mutableData), graph, data,
                                                             mutableData, data.atoms.includeInEachOverlapCheck);
                }
            }
            restoreSidechains(rng, graph, data, mutableData, glycositePreferences,
                              includedGlycanIndices(data, mutableData));
            gmml::log(__LINE__, __FILE__, gmml::INF, "Overlap: " + std::to_string(currentState.overlap.count));
        };

        auto writeOffFile =
            [](const assembly::Graph& graph, const AssemblyData& data, const std::vector<cds::Coordinate>& coordinates,
               const std::vector<bool>& includedMolecules, const std::string& outputDir, const std::string& prefix)
        {
            std::vector<bool> residueIncluded(graph.residueCount, false);
            for (size_t n = 0; n < graph.residueCount; n++)
            {
                residueIncluded[n] = includedMolecules[graph.residueMolecule[n]];
            }
            cds::OffFileData offData = toOffFileData(graph, data, coordinates);
            std::string fileName     = outputDir + "/" + prefix + ".off";
            std::ofstream outFileStream;
            outFileStream.open(fileName.c_str());
            cds::WriteResiduesTogetherToOffFile(outFileStream, graph, offData,
                                                codeUtils::boolsToIndices(residueIncluded), "GLYCOPROTEINBUILDER");
            outFileStream.close();
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
                                const std::vector<cds::Coordinate>& coordinates, const std::vector<int>& atomNumbers,
                                const std::vector<int>& residueNumbers, const std::vector<bool>& includedMolecules,
                                const std::vector<std::array<size_t, 2>>& connectionIndices,
                                const std::string& outputDir, const std::string& prefix)
        {
            cds::PdbFileData pdbData =
                toPdbFileData(graph, data, coordinates, atomNumbers, residueNumbers, headerLines);
            std::vector<std::vector<size_t>> residueIndices =
                codeUtils::boolsToValues(graph.molecules.nodes.constituents, includedMolecules);
            std::vector<std::vector<bool>> TER = residueTER(data, residueIndices);
            codeUtils::createDirectories(outputDir);
            std::string fileName = outputDir + "/" + prefix + ".pdb";
            std::ofstream outFileStream;
            outFileStream.open(fileName.c_str());
            cds::writeAssemblyToPdb(outFileStream, graph, residueIndices, TER, connectionIndices, pdbData);
            outFileStream.close();
        };

        auto printDihedralAnglesAndOverlapOfGlycosites =
            [](const assembly::Graph& graph, const AssemblyData& data, const MutableData& mutableData)
        {
            const std::vector<bool>& included         = glycanIncluded(data, mutableData);
            const std::vector<double>& residueWeights = data.equalResidueWeight;
            const std::vector<bool>& includedAtoms    = data.atoms.includeInEachOverlapCheck;
            size_t glycanCount                        = data.glycans.moleculeId.size();
            for (size_t n = 0; n < glycanCount; n++)
            {
                std::stringstream logss;
                std::string residueID = data.residues.ids[data.glycans.attachmentResidue[n]];
                if (included[n])
                {
                    cds::Overlap selfOverlap = cds::overlapVectorSum(
                        intraGlycanOverlaps(graph, data, mutableData.bounds, residueWeights, includedAtoms, n));
                    cds::Overlap proteinOverlap {0.0, 0.0};
                    for (size_t k : data.indices.proteinMolecules)
                    {
                        proteinOverlap +=
                            cds::overlapVectorSum(moleculeOverlaps(graph, data, mutableData.bounds, residueWeights,
                                                                   includedAtoms, k, data.glycans.moleculeId[n]));
                    }
                    cds::Overlap glycanOverlap {0.0, 0.0};
                    for (size_t k = 0; k < data.glycans.moleculeId.size(); k++)
                    {
                        if (included[k] && k != n)
                        {
                            glycanOverlap += cds::overlapVectorSum(
                                moleculeOverlaps(graph, data, mutableData.bounds, residueWeights, includedAtoms,
                                                 data.glycans.moleculeId[n], data.glycans.moleculeId[k]));
                        }
                    }
                    logss << "Residue ID: " << residueID << ", protein overlap: " << proteinOverlap.count
                          << ", glycan overlap: " << glycanOverlap.count << ", self overlap: " << selfOverlap.count;
                }
                else
                {
                    logss << "Residue ID: " << residueID << ", glycan deleted in order to resolve overlaps";
                }
                gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
            }
        };

        const assembly::Graph& graph    = assembly.graph;
        const AssemblyData& data        = assembly.data;
        const MutableData& initialState = assembly.mutableData;
        SidechainAdjustment sidechainAdjustment =
            settings.allowSidechainAdjustment ? adjustSidechains : noSidechainAdjustment;
        SidechainAdjustment sidechainRestoration =
            settings.allowSidechainAdjustment ? restoreSidechains : noSidechainAdjustment;

        std::vector<cds::ResidueType> nonProteinTypes {cds::ResidueType::Sugar, cds::ResidueType::Derivative,
                                                       cds::ResidueType::Aglycone, cds::ResidueType::Undefined};
        auto isNonProteinResidue = [&data, &nonProteinTypes](size_t n)
        {
            return codeUtils::contains(nonProteinTypes, data.residues.types[n]);
        };

        std::vector<std::array<size_t, 2>> noConnections = {};
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

        auto runInitial = [&](pcg32 rng)
        {
            MutableData mutableData = initialState;
            resolveOverlapsWithWiggler(rng, sidechainAdjustment, sidechainRestoration, graph, data, mutableData, false);
            std::vector<cds::Coordinate> resolvedCoords = getCoordinates(mutableData.bounds.atoms);
            printDihedralAnglesAndOverlapOfGlycosites(graph, data, mutableData);
            writePdbFile(graph, data, resolvedCoords, data.atoms.numbers, data.residues.numbers,
                         mutableData.moleculeIncluded, noConnections, outputDir, "glycoprotein");
            if (settings.writeOffFile)
            {
                writeOffFile(graph, data, resolvedCoords, mutableData.moleculeIncluded, outputDir, "glycoprotein");
            }
            writePdbFile(graph, data, resolvedCoords, data.atoms.serializedNumbers, data.residues.serializedNumbers,
                         mutableData.moleculeIncluded, atomPairsConnectingNonProteinResidues, outputDir,
                         "glycoprotein_serialized");
        };

        auto runIteration = [&](pcg32 rng, size_t count)
        {
            MutableData mutableData = initialState;
            resolveOverlapsWithWiggler(rng, sidechainAdjustment, sidechainRestoration, graph, data, mutableData,
                                       settings.deleteSitesUntilResolved);
            std::vector<cds::Coordinate> coordinates = getCoordinates(mutableData.bounds.atoms);
            bool hasDeleted                          = codeUtils::contains(mutableData.moleculeIncluded, false);
            bool reject                              = false;
            std::string directory                    = outputDir;
            if (settings.rejectExcessiveGlycanOverlaps)
            {
                std::vector<cds::Overlap> atomOverlaps =
                    totalOverlaps(graph, data, mutableData, data.equalResidueWeight,
                                  data.atoms.includeInEachOverlapCheck, OverlapMultiplier {1.0, 1.0, 1.0});
                for (size_t molecule : includedGlycanMoleculeIds(data, mutableData))
                {
                    for (size_t residue : moleculeResidues(graph, molecule))
                    {
                        cds::Overlap overlaps = cds::overlapVectorSum(
                            codeUtils::indicesToValues(atomOverlaps, residueAtoms(graph, residue)));
                        reject = reject || ((overlaps.count * 2.0) >= settings.glycanOverlapRejectionThreshold);
                    }
                }
            }
            directory = hasDeleted ? deletionDir : (reject ? rejectDir : outputDir);
            printDihedralAnglesAndOverlapOfGlycosites(graph, data, mutableData);
            std::stringstream prefix;
            prefix << count << "_glycoprotein";
            std::string prefixStr = prefix.str();
            writePdbFile(graph, data, coordinates, data.atoms.serializedNumbers, data.residues.serializedNumbers,
                         mutableData.moleculeIncluded, atomPairsConnectingNonProteinResidues, directory, prefixStr);
            if (settings.writeOffFile)
            {
                writeOffFile(graph, data, coordinates, mutableData.moleculeIncluded, directory, prefixStr);
            }
            size_t moved = 0;
            for (size_t n : graph.residues.nodes.indices)
            {
                if (mutableData.residueSidechainMoved[n])
                {
                    moved++;
                }
            }
        };

        writePdbFile(graph, data, getCoordinates(data.atoms.initialState), data.atoms.numbers, data.residues.numbers,
                     initialState.moleculeIncluded, noConnections, outputDir, "glycoprotein_initial");

        // order of parallel for loop is undefined, so we generate all seeds up front for determinism
        size_t totalStructures = settings.number3DStructures + 1;
        std::vector<uint64_t> rngSeeds;
        rngSeeds.reserve(totalStructures);
        for (size_t n = 0; n < totalStructures; n++)
        {
            rngSeeds.push_back(codeUtils::generateRandomSeed(seedingRng));
        }

        codeUtils::setOpenMpNumberOfThreads(numThreads);
// clang format doesn't align pragmas
/*    */#pragma omp parallel for
        for (size_t count = 0; count < totalStructures; count++)
        {
            pcg32 rng(rngSeeds[count]);
            if (count == 0)
            {
                runInitial(rng);
            }
            else
            {
                runIteration(rng, count - 1);
            }
        }
    }
} // namespace glycoproteinBuilder
