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
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileWriter.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
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
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/random.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/Graph/types.hpp"
#include "includes/Graph/manipulation.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"
#include "includes/External_Libraries/PCG/pcg_extras.h"

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
                glycosites.emplace_back(glycositeResidue, carb, highestResidueNumber);
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

    void GlycoproteinBuilder::ResolveOverlaps(std::string outputDir)
    {
        double selfWeight           = 1000000.0;
        double proteinWeight        = 1000.0;
        double glycanWeight         = 1.0;
        OverlapWeight overlapWeight = {proteinWeight, glycanWeight, selfWeight};

        uint64_t seed = settings.isDeterministic ? settings.seed : codeUtils::generateRandomSeed();
        pcg32 rng(seed);

        auto randomMetadata = [&rng](GlycamMetadata::DihedralAngleDataVector metadataVector)
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
        auto randomAngle = [&rng, &standardDeviation, &preferenceDeviation](GlycamMetadata::DihedralAngleData metadata)
        {
            double stdCutoff = preferenceDeviation;
            double num       = codeUtils::normalDistributionRandomDoubleWithCutoff(rng, -stdCutoff, stdCutoff);
            auto std         = standardDeviation(metadata);
            return metadata.default_angle + num * (num < 0 ? std.first : std.second);
        };
        auto randomizeShape = [&randomMetadata, &randomAngle, &freezeGlycositeResidueConformation](
                                  const AssemblyGraphs& graphs, const AssemblyData& data, size_t linkageId)
        {
            return randomLinkageShapePreference(graphs, data, linkageId, randomMetadata, randomAngle,
                                                freezeGlycositeResidueConformation);
        };
        auto extractCoordinates = [](const AssemblyData& data)
        {
            const std::vector<cds::Sphere>& bounds = data.atoms.bounds;
            std::vector<cds::Coordinate> result;
            result.reserve(bounds.size());
            for (size_t n = 0; n < bounds.size(); n++)
            {
                result.push_back(bounds[n].center);
            }
            return result;
        };

        auto restoreAllGlycans = [](AssemblyData& data)
        {
            codeUtils::fill(data.glycanData.included, true);
        };

        auto deleteGlycan = [](AssemblyData& data, size_t glycanId)
        {
            data.glycanData.included[glycanId] = false;
        };

        auto resolveOverlapsWithWiggler = [&](const AssemblyGraphs& graphs, AssemblyData& data,
                                              const std::vector<cds::Coordinate>& initialCoordinates,
                                              bool deleteSitesUntilResolved)
        {
            std::vector<std::vector<cds::ResidueLinkageShapePreference>> glycositePreferences;
            for (size_t glycanId = 0; glycanId < graphs.glycans.size(); glycanId++)
            {
                auto preference                       = randomizeShape(graphs, data, glycanId);
                const std::vector<size_t>& linkageIds = graphs.glycans[glycanId].linkages;
                for (size_t k = 0; k < linkageIds.size(); k++)
                {
                    setLinkageShapeToPreference(graphs, data, linkageIds[k], preference[k]);
                }
                updateGlycanBounds(graphs, data, glycanId);
                glycositePreferences.push_back(preference);
            }
            for (size_t glycanId : codeUtils::shuffleVector(rng, codeUtils::indexVector(graphs.glycans)))
            {
                wiggleGlycan(graphs, data, glycanId, searchSettings, overlapWeight, glycositePreferences[glycanId]);
            }
            GlycoproteinState currentState;
            std::vector<size_t> overlapSites =
                determineSitesWithOverlap(codeUtils::indexVector(graphs.glycans), graphs, data);
            for (bool done = false; !done; done = overlapSites.empty() || !deleteSitesUntilResolved)
            {
                cds::Overlap initialOverlap = totalOverlaps(graphs, data, overlapWeight);

                GlycoproteinState initialState = {initialOverlap, overlapSites, glycositePreferences};
                currentState = randomDescent(rng, randomizeShape, searchSettings, settings.persistCycles, overlapWeight,
                                             graphs, data, initialState);
                overlapSites = currentState.overlapSites;
                if (deleteSitesUntilResolved && !overlapSites.empty())
                {
                    size_t indexToRemove = codeUtils::randomIndex(rng, overlapSites);
                    size_t glycan        = overlapSites[indexToRemove];
                    deleteGlycan(data, glycan);
                    overlapSites.erase(overlapSites.begin() + indexToRemove);
                    size_t proteinResidue = graphs.glycans[glycan].attachmentResidue;
                    // restore atoms to initial shape
                    for (size_t n : residueAtoms(graphs, proteinResidue))
                    {
                        data.atoms.bounds[n].center = initialCoordinates[n];
                    }
                    updateResidueBounds(graphs, data, proteinResidue);
                    updateResidueMoleculeBounds(graphs, data, proteinResidue);
                }
            }
            gmml::log(__LINE__, __FILE__, gmml::INF, "Overlap: " + std::to_string(currentState.overlap.count));
            return extractCoordinates(data);
        };

        auto writeOffFile = [](const cds::OffFileData& data, const std::string& prefix)
        {
            std::string fileName = prefix + ".off";
            std::ofstream outFileStream;
            outFileStream.open(fileName.c_str());
            cds::WriteResiduesTogetherToOffFile(outFileStream, data, "GLYCOPROTEINBUILDER");
            outFileStream.close();
        };

        auto writePdbFile = [](const std::vector<std::vector<size_t>>& residueIndices,
                               const std::vector<std::vector<bool>>& residueTER, const cds::PdbFileData& data,
                               const std::vector<std::pair<size_t, size_t>>& connectionIndices,
                               const std::string& prefix)
        {
            std::string fileName = prefix + ".pdb";
            std::ofstream outFileStream;
            outFileStream.open(fileName.c_str());
            for (size_t n = 0; n < residueIndices.size(); n++)
            {
                cds::writeMoleculeToPdb(outFileStream, residueIndices[n], residueTER[n], data);
            }
            cds::writeConectCards(outFileStream, data.atoms.numbers, connectionIndices);
            outFileStream.close();
        };

        auto printDihedralAnglesAndOverlapOfGlycosites = [](const AssemblyGraphs& graphs, const AssemblyData& data)
        {
            const std::vector<GlycanIndices>& glycans = graphs.glycans;
            const std::vector<bool>& included         = data.glycanData.included;
            for (size_t n = 0; n < glycans.size(); n++)
            {
                std::stringstream logss;
                const GlycanIndices& glycan = graphs.glycans[n];
                std::string residueID       = graphs.indices.residues[glycan.attachmentResidue]->getStringId();
                if (included[n])
                {
                    cds::Overlap selfOverlap = intraGlycanOverlaps(graphs, data, n);
                    cds::Overlap proteinOverlap {0.0, 0.0};
                    for (size_t k : graphs.proteinMolecules)
                    {
                        proteinOverlap += moleculeOverlaps(graphs, data, k, glycan.glycanMolecule);
                    }
                    cds::Overlap glycanOverlap {0.0, 0.0};
                    for (size_t k = 0; k < graphs.glycans.size(); k++)
                    {
                        if (included[k] && k != n)
                        {
                            glycanOverlap +=
                                moleculeOverlaps(graphs, data, glycan.glycanMolecule, glycans[k].glycanMolecule);
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

        GlycoproteinAssembly assembly        = toGlycoproteinAssemblyStructs(molecules, glycosites_, overlapWeight);
        AssemblyGraphs& graphs               = assembly.graphs;
        std::vector<cds::Residue*>& residues = graphs.indices.residues;
        std::vector<cds::Atom*>& atoms       = graphs.indices.atoms;
        std::vector<std::vector<size_t>>& moleculeResidues = graphs.molecules.nodes.elements;
        AssemblyData& data                                 = assembly.data;

        std::vector<cds::ResidueType> residueTypes = cds::residueTypes(residues);

        auto includedMolecules = [&](const std::vector<bool>& includedGlycans)
        {
            std::vector<bool> moleculeIncluded(graphs.proteinMolecules.size(), true);
            codeUtils::insertInto(moleculeIncluded, includedGlycans);
            return moleculeIncluded;
        };

        auto residueTER = [&](const std::vector<std::vector<size_t>>& moleculeResidues)
        {
            std::vector<std::vector<bool>> result;
            result.reserve(moleculeResidues.size());
            for (auto& indices : moleculeResidues)
            {
                result.push_back(cds::residueTER(codeUtils::indexValues(residueTypes, indices)));
            }
            return result;
        };

        std::vector<cds::ResidueType> nonProteinTypes {cds::ResidueType::Sugar, cds::ResidueType::Derivative,
                                                       cds::ResidueType::Aglycone, cds::ResidueType::Undefined};
        auto isNonProteinResidue = [&data, &nonProteinTypes](size_t n)
        {
            return codeUtils::contains(nonProteinTypes, data.residues.types[n]);
        };

        std::vector<std::pair<size_t, size_t>> atomPairsConnectingNonProteinResidues;
        const graph::GraphEdges& residueEdges = graphs.residues.edges;
        for (size_t n = 0; n < residueEdges.indices.size(); n++)
        {
            const std::array<size_t, 2>& adj = residueEdges.nodeAdjacencies[n];
            if (isNonProteinResidue(adj[0]) || isNonProteinResidue(adj[1]))
            {
                size_t atomEdgeId                   = residueEdges.indices[n];
                const std::array<size_t, 2> atomAdj = graphs.atoms.edges.nodeAdjacencies[atomEdgeId];
                atomPairsConnectingNonProteinResidues.push_back({atomAdj[0], atomAdj[1]});
            }
        }

        std::vector<std::pair<size_t, size_t>> noConnections = {};

        std::vector<std::vector<bool>> allResidueTER   = residueTER(moleculeResidues);
        std::vector<cds::Coordinate> initalCoordinates = extractCoordinates(data);
        cds::PdbFileData pdbData                       = toPdbFileData(graphs, data);
        writePdbFile(moleculeResidues, allResidueTER, pdbData, noConnections, outputDir + "glycoprotein_initial");
        std::vector<cds::Coordinate> resolvedCoords =
            resolveOverlapsWithWiggler(graphs, data, initalCoordinates, false);
        pdbData.atoms.coordinates = resolvedCoords;
        printDihedralAnglesAndOverlapOfGlycosites(graphs, data);
        writePdbFile(moleculeResidues, allResidueTER, pdbData, noConnections, outputDir + "glycoprotein");
        data.atoms.numbers       = cds::serializedNumberVector(atoms.size());
        data.residues.numbers    = cds::serializedNumberVector(residues.size());
        pdbData.atoms.numbers    = data.atoms.numbers;
        pdbData.residues.numbers = data.residues.numbers;
        {
            cds::OffFileData offData  = toOffFileData(graphs, data);
            offData.atoms.coordinates = resolvedCoords;
            writeOffFile(offData, outputDir + "glycoprotein");
        }
        writePdbFile(moleculeResidues, allResidueTER, pdbData, atomPairsConnectingNonProteinResidues,
                     outputDir + "glycoprotein_serialized");

        for (size_t count = 0; count < settings.number3DStructures; count++)
        {
            restoreAllGlycans(data);
            pdbData.atoms.coordinates =
                resolveOverlapsWithWiggler(graphs, data, initalCoordinates, settings.deleteSitesUntilResolved);
            printDihedralAnglesAndOverlapOfGlycosites(graphs, data);
            std::stringstream prefix;
            prefix << count << "_glycoprotein";
            std::vector<bool> currentMolecules = includedMolecules(data.glycanData.included);
            std::vector<std::vector<size_t>> currentMoleculeResidues =
                codeUtils::maskValues(graphs.molecules.nodes.elements, currentMolecules);
            writePdbFile(currentMoleculeResidues, residueTER(currentMoleculeResidues), pdbData,
                         atomPairsConnectingNonProteinResidues, outputDir + prefix.str());
        }
    }

} // namespace glycoproteinBuilder
