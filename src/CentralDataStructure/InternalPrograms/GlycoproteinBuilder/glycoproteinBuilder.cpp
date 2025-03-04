#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinBuilder.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/cdsInterface.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/writerInterface.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapCount.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/shapeRandomization.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/randomDescent.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileWriter.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileWriter.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageCreation.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
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
    namespace
    {
        Residue* selectResidueFromInput(cds::Assembly* glycoprotein, const std::string userSelection)
        { // Chain_residueNumber_insertionCode* *optional.
            std::vector<std::string> splitUserSelection = codeUtils::split(userSelection, '_');
            if (splitUserSelection.size() < 2)
            {
                throw std::runtime_error(
                    "userSelection (" + userSelection +
                    ") for residue to glycosylate is incorrect format.\nMust be "
                    "chain_residueNumber_insertionCode.\nInsertionCode is optional. Chain can be ? if no "
                    "chain numbers are in input.\nExamples: ?_24_? or ?_24 will use the first residue it "
                    "encounters numbered 24. A_24_B is A chain, residue 24, insertion code B");
            }
            std::string userSelectedChain = "";
            if (splitUserSelection.at(0) != "?")
            {
                userSelectedChain = splitUserSelection.at(0);
            }
            std::string userSelectedResidue = splitUserSelection.at(1);
            if (splitUserSelection.size() == 3)
            {
                userSelectedResidue += splitUserSelection.at(2); // So will be 24A or just 24.
            }
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "We working with " + userSelectedChain + "_" + userSelectedResidue);
            for (auto& residue : glycoprotein->getResidues())
            {
                if (residue->GetType() != cds::ResidueType::Protein)
                {
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              "Glycosite selector skipping non-protein residue: " + residue->getStringId());
                    continue;
                }
                pdb::PdbResidue* pdbResidue = codeUtils::erratic_cast<pdb::PdbResidue*>(residue);
                //            std::cout << pdbResidue->getChainId() << "_";
                //            std::cout << pdbResidue->getNumberAndInsertionCode() << "\n";
                if ((pdbResidue->getChainId() == userSelectedChain) &&
                    (pdbResidue->getNumberAndInsertionCode() == userSelectedResidue))
                {
                    gmml::log(__LINE__, __FILE__, gmml::INF, "Id of selected glycosite: " + pdbResidue->printId());
                    return residue;
                }
            }
            return nullptr;
        }

        std::vector<GlycosylationSite>
        createGlycosites(cds::Assembly* glycoprotein,
                         std::vector<glycoproteinBuilder::GlycositeInput> glycositesInputVector)
        {
            std::vector<GlycosylationSite> glycosites;
            for (auto& glycositeInput : glycositesInputVector)
            {
                gmml::log(__LINE__, __FILE__, gmml::INF,
                          "Creating glycosite on residue " + glycositeInput.proteinResidueId + " with glycan " +
                              glycositeInput.glycanInput);
                Residue* glycositeResidue = selectResidueFromInput(glycoprotein, glycositeInput.proteinResidueId);
                if (glycositeResidue == nullptr)
                {
                    throw std::runtime_error("Error: Did not find a residue with id matching " +
                                             glycositeInput.proteinResidueId + "\n");
                }
                Carbohydrate* carb = codeUtils::erratic_cast<Carbohydrate*>(
                    glycoprotein->addMolecule(std::make_unique<Carbohydrate>(glycositeInput.glycanInput)));
                gmml::log(__LINE__, __FILE__, gmml::INF,
                          "About to emplace_back to glycosites with: " + glycositeInput.proteinResidueId +
                              " and glycan " + glycositeInput.glycanInput);
                unsigned int highestResidueNumber =
                    cdsSelections::findHighestResidueNumber(glycoprotein->getResidues());
                glycosites.emplace_back(glycositeResidue, carb, glycositeInput, highestResidueNumber);
                for (auto& linkage : carb->GetGlycosidicLinkages())
                {
                    cds::determineResiduesForOverlapCheck(linkage); // Now that the protein residue is attached.
                }
                //	    std::cout << "Done with glycan" << std::endl;
                gmml::log(__LINE__, __FILE__, gmml::INF,
                          "Completed creating glycosite on residue " + glycositeInput.proteinResidueId +
                              " with glycan " + glycositeInput.glycanInput);
            }
            return glycosites;
        }
    } // namespace

    GlycoproteinBuilder::GlycoproteinBuilder(glycoproteinBuilder::GlycoproteinBuilderInputs inputStruct,
                                             pdb::PreprocessorOptions preprocessingOptions)
        : settings(inputStruct)
    {
        try
        {
            pdbFile = pdb::PdbFile(inputStruct.substrateFileName);
            if (!inputStruct.skipMDPrep)
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, "Performing MDPrep aka preprocessing.");
                pdbFile.PreProcess(preprocessingOptions);
            }
            glycoprotein_                                = &pdbFile.mutableAssemblies().front();
            std::vector<cds::Residue*> gpInitialResidues = glycoprotein_->getResidues();
            cds::setIntraConnectivity(gpInitialResidues);
            gmml::log(__LINE__, __FILE__, gmml::INF, "Attaching Glycans To Glycosites.");
            proteinResidues_ = glycoprotein_->getResidues();
            glycosites_      = createGlycosites(glycoprotein_, inputStruct.glycositesInputVector);
            cds::setInterConnectivity(gpInitialResidues); // do the inter here, so that the whole protein isn't included
                                                          // as overlap residues in the glycan linkages.
        }
        catch (const std::string& errorMessage)
        {
            gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
            throw std::runtime_error(errorMessage);
        }
        gmml::log(__LINE__, __FILE__, gmml::INF, "Initialization of Glycoprotein builder complete!");
        return;
    }

    void GlycoproteinBuilder::ResolveOverlaps(const std::string& outputDir, const std::vector<std::string>& headerLines,
                                              int numThreads)
    {
        double selfWeight           = 1000000.0;
        double proteinWeight        = 1000.0;
        double glycanWeight         = 1.0;
        OverlapWeight overlapWeight = {proteinWeight, glycanWeight, selfWeight};

        uint64_t seed = settings.isDeterministic ? settings.seed : codeUtils::generateRandomSeed();
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
            return randomLinkageShapePreference(rng, data, mutableData, linkageId, randomMetadata, randomAngle,
                                                freezeGlycositeResidueConformation);
        };

        auto atomCoordinates = [](const MutableData& data)
        {
            const std::vector<cds::Sphere>& bounds = data.atomBounds;
            std::vector<cds::Coordinate> result;
            result.reserve(bounds.size());
            for (size_t n = 0; n < bounds.size(); n++)
            {
                result.push_back(bounds[n].center);
            }
            return result;
        };

        auto deleteGlycan = [](MutableData& mutableData, size_t glycanId)
        {
            mutableData.glycanIncluded[glycanId] = false;
        };

        auto resolveOverlapsWithWiggler =
            [&](pcg32& rng, const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                const std::vector<cds::Coordinate>& initialCoordinates, bool deleteSitesUntilResolved)
        {
            std::vector<std::vector<cds::ResidueLinkageShapePreference>> glycositePreferences;
            for (size_t glycanId = 0; glycanId < data.indices.glycans.size(); glycanId++)
            {
                auto preference                       = randomizeShape(rng, data, mutableData, glycanId);
                const std::vector<size_t>& linkageIds = data.indices.glycans[glycanId].linkages;
                for (size_t k = 0; k < linkageIds.size(); k++)
                {
                    setLinkageShapeToPreference(graph, data, mutableData, linkageIds[k], preference[k]);
                }
                updateGlycanBounds(graph, data, mutableData, glycanId);
                glycositePreferences.push_back(preference);
            }
            for (size_t glycanId : codeUtils::shuffleVector(rng, codeUtils::indexVector(data.indices.glycans)))
            {
                wiggleGlycan(graph, data, mutableData, glycanId, searchSettings, overlapWeight,
                             glycositePreferences[glycanId]);
            }
            GlycoproteinState currentState;
            std::vector<size_t> overlapSites =
                determineSitesWithOverlap(codeUtils::indexVector(data.indices.glycans), graph, data, mutableData);
            for (bool done = false; !done; done = overlapSites.empty() || !deleteSitesUntilResolved)
            {
                cds::Overlap initialOverlap = totalOverlaps(graph, data, mutableData, overlapWeight);

                GlycoproteinState initialState = {initialOverlap, overlapSites, glycositePreferences};
                currentState = randomDescent(rng, randomizeShape, searchSettings, settings.persistCycles, overlapWeight,
                                             graph, data, mutableData, initialState);
                overlapSites = currentState.overlapSites;
                if (deleteSitesUntilResolved && !overlapSites.empty())
                {
                    size_t indexToRemove = codeUtils::randomIndex(rng, overlapSites);
                    size_t glycan        = overlapSites[indexToRemove];
                    deleteGlycan(mutableData, glycan);
                    size_t proteinResidue = data.indices.glycans[glycan].attachmentResidue;
                    // restore atoms to initial shape
                    for (size_t n : residueAtoms(graph, proteinResidue))
                    {
                        mutableData.atomBounds[n].center = initialCoordinates[n];
                    }
                    updateResidueBounds(graph, mutableData, proteinResidue);
                    updateResidueMoleculeBounds(graph, mutableData, proteinResidue);
                    overlapSites = determineSitesWithOverlap(codeUtils::indexVector(data.indices.glycans), graph, data,
                                                             mutableData);
                }
            }
            gmml::log(__LINE__, __FILE__, gmml::INF, "Overlap: " + std::to_string(currentState.overlap.count));
            return atomCoordinates(mutableData);
        };

        auto writeOffFile = [&outputDir](const assembly::Graph& graph, const AssemblyData& data,
                                         const std::vector<cds::Coordinate>& coordinates, const std::string& prefix)
        {
            cds::OffFileData offData = toOffFileData(graph, data, coordinates);
            std::string fileName     = outputDir + "/" + prefix + ".off";
            std::ofstream outFileStream;
            outFileStream.open(fileName.c_str());
            cds::WriteResiduesTogetherToOffFile(outFileStream, graph, offData, "GLYCOPROTEINBUILDER");
            outFileStream.close();
        };

        auto writePdbFile = [&outputDir, &headerLines](
                                const assembly::Graph& graph, const AssemblyData& data,
                                const std::vector<cds::Coordinate>& coordinates, const std::vector<int>& atomNumbers,
                                const std::vector<int>& residueNumbers,
                                const std::vector<std::vector<size_t>>& residueIndices,
                                const std::vector<std::vector<bool>>& residueTER,
                                const std::vector<std::array<size_t, 2>>& connectionIndices, const std::string& prefix)
        {
            cds::PdbFileData pdbData =
                toPdbFileData(graph, data, coordinates, atomNumbers, residueNumbers, headerLines);
            std::string fileName = outputDir + "/" + prefix + ".pdb";
            std::ofstream outFileStream;
            outFileStream.open(fileName.c_str());
            cds::writeAssemblyToPdb(outFileStream, graph, residueIndices, residueTER, connectionIndices, pdbData);
            cds::theEnd(outFileStream);
            outFileStream.close();
        };

        auto printDihedralAnglesAndOverlapOfGlycosites =
            [](const assembly::Graph& graph, const AssemblyData& data, const MutableData& mutableData)
        {
            const std::vector<GlycanIndices>& glycans = data.indices.glycans;
            const std::vector<bool>& included         = mutableData.glycanIncluded;
            for (size_t n = 0; n < glycans.size(); n++)
            {
                std::stringstream logss;
                const GlycanIndices& glycan = data.indices.glycans[n];
                std::string residueID       = data.residues.ids[glycan.attachmentResidue];
                if (included[n])
                {
                    cds::Overlap selfOverlap = intraGlycanOverlaps(graph, data, mutableData, n);
                    cds::Overlap proteinOverlap {0.0, 0.0};
                    for (size_t k : data.indices.proteinMolecules)
                    {
                        proteinOverlap += moleculeOverlaps(graph, data, mutableData, k, glycan.glycanMolecule);
                    }
                    cds::Overlap glycanOverlap {0.0, 0.0};
                    for (size_t k = 0; k < data.indices.glycans.size(); k++)
                    {
                        if (included[k] && k != n)
                        {
                            glycanOverlap += moleculeOverlaps(graph, data, mutableData, glycan.glycanMolecule,
                                                              glycans[k].glycanMolecule);
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

        std::vector<cds::Molecule*> molecules = getGlycoprotein()->getMolecules();

        const GlycoproteinAssembly assembly = toGlycoproteinAssemblyStructs(molecules, glycosites_, overlapWeight);
        const assembly::Graph& graph        = assembly.graph;
        const std::vector<std::vector<size_t>>& moleculeResidues = graph.molecules.nodes.elements;
        const AssemblyData& data                                 = assembly.data;
        const MutableData& initialState                          = assembly.mutableData;
        std::vector<cds::Coordinate> initialCoordinates          = atomCoordinates(initialState);

        auto includedMolecules = [&](const std::vector<bool>& includedGlycans)
        {
            std::vector<bool> proteinIncluded(data.indices.proteinMolecules.size(), true);
            return codeUtils::vectorAppend(proteinIncluded, includedGlycans);
        };

        auto residueTER = [&](const std::vector<std::vector<size_t>>& moleculeResidues)
        {
            std::vector<std::vector<bool>> result;
            result.reserve(moleculeResidues.size());
            for (auto& indices : moleculeResidues)
            {
                result.push_back(cds::residueTER(codeUtils::indicesToValues(data.residues.types, indices)));
            }
            return result;
        };

        std::vector<cds::ResidueType> nonProteinTypes {cds::ResidueType::Sugar, cds::ResidueType::Derivative,
                                                       cds::ResidueType::Aglycone, cds::ResidueType::Undefined};
        auto isNonProteinResidue = [&data, &nonProteinTypes](size_t n)
        {
            return codeUtils::contains(nonProteinTypes, data.residues.types[n]);
        };

        std::vector<std::array<size_t, 2>> noConnections = {};
        std::vector<std::vector<bool>> allResidueTER     = residueTER(moleculeResidues);
        std::vector<std::array<size_t, 2>> atomPairsConnectingNonProteinResidues;
        const graph::GraphEdges& residueEdges = graph.residues.edges;
        for (size_t n = 0; n < residueEdges.indices.size(); n++)
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
            std::vector<cds::Coordinate> resolvedCoords =
                resolveOverlapsWithWiggler(rng, graph, data, mutableData, initialCoordinates, false);
            printDihedralAnglesAndOverlapOfGlycosites(graph, data, mutableData);
            writePdbFile(graph, data, resolvedCoords, data.atoms.numbers, data.residues.numbers, moleculeResidues,
                         allResidueTER, noConnections, "glycoprotein");
            writeOffFile(graph, data, resolvedCoords, "glycoprotein");
            writePdbFile(graph, data, resolvedCoords, data.atoms.serializedNumbers, data.residues.serializedNumbers,
                         moleculeResidues, allResidueTER, atomPairsConnectingNonProteinResidues,
                         "glycoprotein_serialized");
        };

        auto runIteration = [&](pcg32 rng, size_t count)
        {
            MutableData mutableData                  = initialState;
            std::vector<cds::Coordinate> coordinates = resolveOverlapsWithWiggler(
                rng, graph, data, mutableData, initialCoordinates, settings.deleteSitesUntilResolved);
            printDihedralAnglesAndOverlapOfGlycosites(graph, data, mutableData);
            std::stringstream prefix;
            prefix << count << "_glycoprotein";
            std::vector<bool> currentMolecules = includedMolecules(mutableData.glycanIncluded);
            std::vector<std::vector<size_t>> currentMoleculeResidues =
                codeUtils::boolsToValues(graph.molecules.nodes.elements, currentMolecules);
            writePdbFile(graph, data, coordinates, data.atoms.serializedNumbers, data.residues.serializedNumbers,
                         currentMoleculeResidues, residueTER(currentMoleculeResidues),
                         atomPairsConnectingNonProteinResidues, prefix.str());
        };

        writePdbFile(graph, data, initialCoordinates, data.atoms.numbers, data.residues.numbers, moleculeResidues,
                     allResidueTER, noConnections, "glycoprotein_initial");

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
