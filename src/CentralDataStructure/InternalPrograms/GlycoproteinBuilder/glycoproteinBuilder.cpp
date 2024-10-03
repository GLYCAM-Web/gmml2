#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinBuilder.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinOverlapResolution.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
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
#include "includes/CentralDataStructure/Geometry/types.hpp"
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

struct GlycositeData
{
    std::vector<std::vector<cds::ResidueLinkage>> glycosidicLinkages;
    std::vector<std::vector<Residue*>> glycositeResidues;
};

namespace
{
    GlycositeData toGlycositeData(std::vector<GlycosylationSite>& glycosites)
    {
        std::vector<std::vector<cds::ResidueLinkage>> glycosidicLinkages;
        std::vector<std::vector<Residue*>> glycositeResidues;
        for (auto& a : glycosites)
        {
            auto glycan = a.GetGlycan();
            glycosidicLinkages.push_back(glycan->GetGlycosidicLinkages());
            glycositeResidues.push_back(glycan->getResidues());
        }
        return GlycositeData {glycosidicLinkages, glycositeResidues};
    }

    std::vector<OverlapResidues> toOverlapResidues(const std::vector<cds::Residue*>& proteinResidues,
                                                   std::vector<GlycosylationSite>& glycosites,
                                                   const std::vector<std::vector<Residue*>>& glycositeResidues)
    {
        std::vector<OverlapResidues> result;
        for (size_t n = 0; n < glycosites.size(); n++)
        {
            std::vector<Residue*> protein = proteinResidues;
            auto glycosite                = glycosites[n].GetResidue();
            protein.erase(std::remove(protein.begin(), protein.end(), glycosite), protein.end());
            protein.insert(protein.begin(), glycosite);

            std::vector<Residue*> glycan;
            for (auto& otherSite : codeUtils::withoutNth(n, glycositeResidues))
            {
                codeUtils::insertInto(glycan, otherSite);
            }

            result.push_back({protein, glycan});
        }
        return result;
    }

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
        gmml::log(__LINE__, __FILE__, gmml::INF, "We working with " + userSelectedChain + "_" + userSelectedResidue);
        for (auto& residue : glycoprotein->getResidues())
        {
            pdb::PdbResidue* pdbResidue = codeUtils::erratic_cast<pdb::PdbResidue*>(residue);
            // std::cout << pdbResidue->getChainId() << "_";
            //        std::cout << pdbResidue->getNumberAndInsertionCode() << "\n";
            if ((pdbResidue->getChainId() == userSelectedChain) &&
                (pdbResidue->getNumberAndInsertionCode() == userSelectedResidue))
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, "Id of selected glycosite: " + pdbResidue->printId());
                return residue;
            }
        }
        return nullptr;
    }

    std::vector<GlycosylationSite> createGlycosites(cds::Assembly* glycoprotein,
                                                    std::vector<glycoprotein::GlycositeInput> glycositesInputVector)
    {
        std::vector<GlycosylationSite> glycosites;
        for (auto& glycositeInput : glycositesInputVector)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Creating glycosite on residue " + glycositeInput.proteinResidueId + " with glycan " +
                          glycositeInput.glycanInput);
            Carbohydrate* carb = codeUtils::erratic_cast<Carbohydrate*>(
                glycoprotein->addMolecule(std::make_unique<Carbohydrate>(glycositeInput.glycanInput)));
            Residue* glycositeResidue = selectResidueFromInput(glycoprotein, glycositeInput.proteinResidueId);
            if (glycositeResidue == nullptr)
            {
                throw std::runtime_error("Did not find a residue with id matching " + glycositeInput.proteinResidueId);
            }
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "About to emplace_back to glycosites with: " + glycositeInput.proteinResidueId + " and glycan " +
                          glycositeInput.glycanInput);
            unsigned int highestResidueNumber = cdsSelections::findHighestResidueNumber(glycoprotein->getResidues());
            glycosites.emplace_back(glycositeResidue, carb, highestResidueNumber);
            for (auto& linkage : carb->GetGlycosidicLinkages())
            {
                cds::determineResiduesForOverlapCheck(linkage); // Now that the protein residue is attached.
            }
            //	    std::cout << "Done with glycan" << std::endl;
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Completed creating glycosite on residue " + glycositeInput.proteinResidueId + " with glycan " +
                          glycositeInput.glycanInput);
        }
        return glycosites;
    }
} // namespace

GlycoproteinBuilder::GlycoproteinBuilder(glycoprotein::GlycoproteinBuilderInputs inputStruct,
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
        cds::setInterConnectivity(gpInitialResidues); // do the inter here, so that the whole protein isn't included as
                                                      // overlap residues in the glycan linkages.
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
    double selfWeight    = 1000000.0;
    double proteinWeight = 1000.0;
    double glycanWeight  = 1.0;
    auto overlapWeight   = OverlapWeight {proteinWeight, glycanWeight, selfWeight};

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
    auto standardDeviation = [&preferenceDeviation, &searchDeviation](const GlycamMetadata::DihedralAngleData& metadata)
    {
        auto angle = metadata.angle_deviation;
        if (std::holds_alternative<GlycamMetadata::AngleLimit>(angle))
        {
            auto dev         = std::get<GlycamMetadata::AngleLimit>(angle);
            double max_std   = preferenceDeviation + searchDeviation;
            double lower_std = dev.lowerDeviationLimit / max_std;
            double upper_std = dev.upperDeviationLimit / max_std;
            return std::pair<double, double> {lower_std, upper_std};
        }
        else if (std::holds_alternative<GlycamMetadata::AngleStd>(angle))
        {
            auto dev = std::get<GlycamMetadata::AngleStd>(angle);
            return std::pair<double, double> {dev.lowerDeviationStd, dev.upperDeviationStd};
        }
        else
        {
            throw std::runtime_error("unknown angle deviation type");
        }
    };
    auto searchAngles = [&standardDeviation, &searchIncrement](const GlycamMetadata::DihedralAngleData& metadata,
                                                               double preference, double deviation)
    {
        auto std = standardDeviation(metadata);
        return cds::evenlySpacedAngles(preference, deviation * std.first, deviation * std.second, searchIncrement);
    };
    auto searchSettings = cds::AngleSearchSettings {searchDeviation, searchAngles};
    auto randomAngle    = [&rng, &standardDeviation, &preferenceDeviation](GlycamMetadata::DihedralAngleData metadata)
    {
        double stdCutoff = preferenceDeviation;
        double num       = codeUtils::normalDistributionRandomDoubleWithCutoff(rng, -stdCutoff, stdCutoff);
        auto std         = standardDeviation(metadata);
        return metadata.default_angle + num * (num < 0 ? std.first : std.second);
    };
    auto glycositeResidueConformationFrozen =
        [](std::vector<cds::ResidueLinkageShapePreference> preference, const std::vector<cds::ResidueLinkage>& linkages)
    {
        if (std::holds_alternative<cds::ConformerShapePreference>(preference[0]))
        {
            auto& firstLinkage = linkages[0];
            auto& pref         = std::get<cds::ConformerShapePreference>(preference[0]);
            auto shape         = cds::currentShape(firstLinkage.rotatableDihedrals, firstLinkage.dihedralMetadata);
            for (size_t n = 0; n < firstLinkage.rotatableDihedrals.size(); n++)
            {
                auto& name = firstLinkage.dihedralMetadata[n][0].dihedral_angle_name_;
                if ((name == "Chi1") || (name == "Chi2"))
                {
                    pref.isFrozen[n]         = true;
                    size_t metadata          = shape[n].metadataIndex;
                    pref.angles[n][metadata] = shape[n].value;
                    pref.metadataOrder       = {metadata};
                }
            }
        }
        return preference;
    };
    auto randomizeShape = [&randomMetadata, &randomAngle, &freezeGlycositeResidueConformation,
                           &glycositeResidueConformationFrozen](const std::vector<cds::ResidueLinkage>& linkages)
    {
        auto preference = cds::linkageShapePreference(randomMetadata, randomAngle, linkages);
        return freezeGlycositeResidueConformation ? glycositeResidueConformationFrozen(preference, linkages)
                                                  : preference;
    };
    auto wiggleGlycanFunc = [&searchSettings](const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId,
                                              OverlapWeight weight,
                                              const std::vector<cds::ResidueLinkageShapePreference>& preferences)
    {
        return wiggleGlycan(graphs, data, glycanId, searchSettings, weight, preferences);
    };

    auto resolveOverlapsWithWiggler = [&](const AssemblyGraphs& assemblyGraphs, AssemblyData& assemblyData,
                                          GlycositeData& data, const std::vector<cds::Residue*>& proteinResidues,
                                          std::vector<GlycosylationSite>& glycosites)
    {
        auto& glycosidicLinkages = data.glycosidicLinkages;
        auto overlapResidues     = toOverlapResidues(proteinResidues, glycosites, data.glycositeResidues);

        std::vector<std::vector<cds::ResidueLinkageShapePreference>> glycositePreferences;
        std::vector<std::vector<std::vector<cds::AngleWithMetadata>>> glycositeShape;
        for (size_t n = 0; n < glycosidicLinkages.size(); n++)
        {
            auto& linkages                        = glycosidicLinkages[n];
            auto preference                       = randomizeShape(linkages);
            const std::vector<size_t>& linkageIds = assemblyGraphs.glycans[n].linkages;
            for (size_t k = 0; k < linkageIds.size(); k++)
            {
                setLinkageShapeToPreference(assemblyGraphs, assemblyData, linkageIds[k], preference[k]);
            }
            updateGlycanBounds(assemblyGraphs, assemblyData, n);
            glycositePreferences.push_back(preference);
            glycositeShape.push_back(cds::currentShape(linkages));
        }
        for (size_t glycanId : codeUtils::shuffleVector(rng, codeUtils::indexVector(glycosidicLinkages)))
        {
            glycositeShape[glycanId] =
                wiggleGlycanFunc(assemblyGraphs, assemblyData, glycanId, overlapWeight, glycositePreferences[glycanId]);
        }
        cds::Overlap initialOverlap = totalOverlaps(overlapWeight, assemblyGraphs, assemblyData);
        auto overlapSites =
            determineSitesWithOverlap(codeUtils::indexVector(assemblyGraphs.glycans), assemblyGraphs, assemblyData);
        auto initialState = GlycoproteinState {initialOverlap, overlapSites, glycositePreferences, glycositeShape};
        GlycoproteinState currentState =
            randomDescent(rng, randomizeShape, wiggleGlycanFunc, settings.persistCycles, overlapWeight, assemblyGraphs,
                          assemblyData, glycosidicLinkages, initialState);
        for (size_t n = 0; n < assemblyGraphs.indices.atoms.size(); n++)
        {
            assemblyData.atoms.coordinateReferences[n].set(assemblyData.atoms.bounds[n].center);
        }
        gmml::log(__LINE__, __FILE__, gmml::INF, "Overlap: " + std::to_string(currentState.overlap.count));
    };

    auto writeOffFile = [](const cds::OffWriterData& data, const std::string& prefix)
    {
        std::string fileName = prefix + ".off";
        std::ofstream outFileStream;
        outFileStream.open(fileName.c_str());
        cds::WriteResiduesTogetherToOffFile(outFileStream, data, "GLYCOPROTEINBUILDER");
        outFileStream.close();
    };

    auto writePdbFile = [](const std::vector<std::vector<size_t>>& residueIndices,
                           const std::vector<std::vector<bool>>& residueTER, const cds::PdbWriterData& data,
                           const std::vector<std::pair<size_t, size_t>>& connectionIndices, const std::string& prefix)
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
        std::stringstream logss;
        const std::vector<GlycanIndices>& glycans = graphs.glycans;
        for (size_t n = 0; n < glycans.size(); n++)
        {
            const GlycanIndices& glycan = graphs.glycans[n];
            cds::Overlap selfOverlap    = intraGlycanOverlaps(graphs, data, n);
            cds::Overlap proteinOverlap {0.0, 0.0};
            for (size_t n : graphs.proteinMolecules)
            {
                proteinOverlap += moleculeOverlaps(graphs, data, n, glycan.glycanMolecule);
            }
            cds::Overlap glycanOverlap {0.0, 0.0};
            for (size_t k = 0; k < graphs.glycans.size(); k++)
            {
                if (k != n)
                {
                    glycanOverlap += moleculeOverlaps(graphs, data, glycan.glycanMolecule, glycans[k].glycanMolecule);
                }
            }
            logss << "Residue ID: " << graphs.indices.residues[glycan.attachmentResidue]->getStringId()
                  << ", protein overlap: " << proteinOverlap.count << ", glycan overlap: " << glycanOverlap.count
                  << ", self overlap: " << selfOverlap.count;
            gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
        }
    };

    std::vector<cds::Molecule*> molecules         = getGlycoprotein()->getMolecules();
    cds::GraphIndexData graphIndices              = cds::toIndexData(molecules);
    graph::Database atomGraphData                 = cds::createGraphData(graphIndices);
    graph::Graph atomGraph                        = graph::identity(atomGraphData);
    graph::Graph residueGraph                     = graph::quotient(atomGraphData, graphIndices.atomResidue);
    std::vector<std::vector<size_t>>& atomIndices = residueGraph.nodes.elements;
    graph::Graph moleculeGraph = graph::quotient(graph::asData(residueGraph), graphIndices.residueMolecule);
    std::vector<std::vector<size_t>>& residueIndices = moleculeGraph.nodes.elements;
    std::vector<cds::Atom*>& atoms                   = graphIndices.atoms;
    std::vector<cds::Residue*>& residues             = graphIndices.residues;

    std::vector<cds::Coordinate> atomCoordinates = cds::atomCoordinates(atoms);
    std::vector<double> atomRadii                = cds::atomRadii(atoms);
    std::vector<cds::Sphere> atomBoundingSpheres;
    atomBoundingSpheres.reserve(atoms.size());
    for (size_t n = 0; n < atoms.size(); n++)
    {
        atomBoundingSpheres.push_back({atomRadii[n], atomCoordinates[n]});
    }
    AtomData atomData {cds::atomCoordinateReferences(atoms), atomCoordinates, atomBoundingSpheres};
    auto boundingSpheresOf =
        [](const std::vector<cds::Sphere>& spheres, const std::vector<std::vector<size_t>>& indexVector)
    {
        std::vector<cds::Sphere> result;
        result.reserve(indexVector.size());
        for (auto& indices : indexVector)
        {
            result.push_back(cds::boundingSphere(codeUtils::indexValues(spheres, indices)));
        }
        return result;
    };

    auto indexOfResidues = [&residues](const std::vector<Residue*>& toFind)
    {
        std::vector<size_t> result;
        result.reserve(toFind.size());
        for (auto& res : toFind)
        {
            result.push_back(codeUtils::indexOf(residues, res));
        }
        return result;
    };

    auto indexOfAtoms = [&atoms](const std::vector<Atom*>& toFind)
    {
        std::vector<size_t> result;
        result.reserve(toFind.size());
        for (auto& res : toFind)
        {
            result.push_back(codeUtils::indexOf(atoms, res));
        }
        return result;
    };

    auto closelyBondedAtoms = [&](size_t residueId, size_t atomId)
    {
        if (graphIndices.atomResidue[atomId] != residueId)
        {
            throw std::runtime_error("bork");
        }
        const std::vector<size_t>& elements = residueGraph.nodes.elements[residueId];
        std::vector<bool> result(elements.size(), false);
        for (size_t n = 0; n < elements.size(); n++)
        {
            const std::vector<size_t>& adjacencies = atomGraph.nodes.nodeAdjacencies[atomId];
            size_t id                              = elements[n];
            result[n]                              = (id == atomId) || codeUtils::contains(adjacencies, id);
        }
        return result;
    };

    std::vector<MoleculeType> moleculeTypes(graphIndices.molecules.size(), MoleculeType::protein);

    auto glycositeData = toGlycositeData(glycosites_);
    std::vector<size_t> rotatableDihedralCurrentMetadataIndex;
    std::vector<RotatableDihedralIndices> rotatableDihedralIndices;
    std::vector<ResidueLinkageIndices> residueLinkages;
    std::vector<GlycamMetadata::RotamerType> linkageRotamerTypes;
    std::vector<std::vector<GlycamMetadata::DihedralAngleDataVector>> linkageMetadata;
    std::vector<std::vector<cds::BondedResidueOverlapInput>> linkageOverlapBonds;
    std::vector<bool> linkageBranching;
    std::vector<GlycanIndices> glycositeIndices;
    for (size_t n = 0; n < glycosites_.size(); n++)
    {
        size_t site = codeUtils::indexOf(residues, glycosites_[n].GetResidue());
        size_t moleculeIndex =
            codeUtils::indexOf(molecules, codeUtils::erratic_cast<cds::Molecule*>(glycosites_[n].GetGlycan()));
        auto& linkages                 = glycositeData.glycosidicLinkages[n];
        std::vector<size_t> linkageIds = codeUtils::indexVectorWithOffset(residueLinkages.size(), linkages);
        for (auto& linkage : linkages)
        {
            auto& linkageDihedrals = linkage.rotatableDihedrals;
            std::vector<size_t> dihedralIndices =
                codeUtils::indexVectorWithOffset(rotatableDihedralIndices.size(), linkageDihedrals);
            for (auto& dihedral : linkageDihedrals)
            {
                rotatableDihedralCurrentMetadataIndex.push_back(0);
                std::array<size_t, 4> dihedralAtoms;
                for (size_t n = 0; n < 4; n++)
                {
                    dihedralAtoms[n] = codeUtils::indexOf(atoms, dihedral.atoms[n]);
                }
                rotatableDihedralIndices.push_back({dihedralAtoms, indexOfAtoms(dihedral.movingAtoms)});
            }
            auto onlyThisMolecule = [&](const std::vector<size_t>& indices)
            {
                std::vector<size_t> result;
                result.reserve(indices.size());
                for (size_t index : indices)
                {
                    if (graphIndices.residueMolecule[index] == moleculeIndex)
                    {
                        result.push_back(index);
                    }
                }
                return result;
            };
            std::vector<size_t> nonReducing = onlyThisMolecule(indexOfResidues(linkage.nonReducingOverlapResidues));
            std::vector<size_t> reducing    = onlyThisMolecule(indexOfResidues(linkage.reducingOverlapResidues));
            size_t firstResidue             = codeUtils::indexOf(residues, linkage.link.residues.first);
            size_t secondResidue            = codeUtils::indexOf(residues, linkage.link.residues.second);
            auto& adjacencies               = residueGraph.nodes.nodeAdjacencies[firstResidue];
            size_t edgeN                    = codeUtils::indexOf(adjacencies, secondResidue);
            if (edgeN >= adjacencies.size())
            {
                throw std::runtime_error("no residue adjacency");
            }
            size_t edgeId                    = residueGraph.nodes.edgeAdjacencies[firstResidue][edgeN];
            std::array<size_t, 2> residueIds = residueGraph.edges.nodeAdjacencies[edgeId];
            std::array<size_t, 2> atomIds    = atomGraph.edges.nodeAdjacencies[residueGraph.edges.indices[edgeId]];
            bool direction                   = graphIndices.atomResidue[atomIds[0]] == residueIds[1];
            std::array<std::vector<bool>, 2> bondedAtoms = {closelyBondedAtoms(residueIds[0], atomIds[direction]),
                                                            closelyBondedAtoms(residueIds[1], atomIds[!direction])};

            residueLinkages.push_back({edgeId, dihedralIndices, bondedAtoms, nonReducing, reducing});
            linkageRotamerTypes.push_back(linkage.rotamerType);
            linkageMetadata.push_back(linkage.dihedralMetadata);
            linkageOverlapBonds.push_back({
                {residueIds, bondedAtoms}
            });
            linkageBranching.push_back(linkage.rotatableDihedrals[0].isBranchingLinkage);
        }
        moleculeTypes[moleculeIndex] = MoleculeType::glycan;
        glycositeIndices.push_back({site, moleculeIndex, linkageIds});
    }
    std::vector<double> residueOverlapWeight;
    for (size_t molecule : graphIndices.residueMolecule)
    {
        residueOverlapWeight.push_back(moleculeTypes[molecule] == MoleculeType::protein ? proteinWeight : glycanWeight);
    }

    std::vector<cds::Sphere> residueBoundingSpheres =
        boundingSpheresOf(atomBoundingSpheres, residueGraph.nodes.elements);
    ResidueData residueData {residueOverlapWeight, residueBoundingSpheres};
    std::vector<cds::Sphere> moleculeBoundingSpheres =
        boundingSpheresOf(residueBoundingSpheres, moleculeGraph.nodes.elements);
    MoleculeData moleculeData {moleculeTypes, moleculeBoundingSpheres};
    RotatableDihedralData RotatableDihedralData {rotatableDihedralCurrentMetadataIndex};
    ResidueLinkagedata residueLinkageData {linkageRotamerTypes, linkageMetadata, linkageOverlapBonds, linkageBranching};
    AssemblyData assemblyData {atomData, residueData, moleculeData, RotatableDihedralData, residueLinkageData};

    std::vector<size_t> proteinMolecules;
    for (size_t n = 0; n < moleculeTypes.size(); n++)
    {
        if (moleculeTypes[n] == MoleculeType::protein)
        {
            proteinMolecules.push_back(n);
        }
    }

    AssemblyGraphs assemblyGraphs {graphIndices,    atomGraph,        residueGraph,
                                   moleculeGraph,   proteinMolecules, rotatableDihedralIndices,
                                   residueLinkages, glycositeIndices};

    std::vector<cds::ResidueType> residueTypes = cds::residueTypes(residues);
    std::vector<std::vector<bool>> residueTER;

    for (auto& indices : residueIndices)
    {
        residueTER.push_back(cds::residueTER(codeUtils::indexValues(residueTypes, indices)));
    }

    std::vector<std::string> recordNames(atoms.size(), "ATOM");
    std::vector<std::string> chainIds(residues.size(), "");
    std::vector<std::string> insertionCodes(residues.size(), "");
    cds::ResiduePdbData residuePdbData(atomIndices, cds::residueNumbers(residues), cds::residueNames(residues),
                                       chainIds, insertionCodes);
    cds::AtomPdbData atomPdbData(atoms, recordNames);
    cds::PdbWriterData writerData {residuePdbData, atomPdbData};

    auto pdbResidues =
        cdsSelections::selectResiduesByType(residues, {cds::ResidueType::Sugar, cds::ResidueType::Derivative,
                                                       cds::ResidueType::Aglycone, cds::ResidueType::Undefined});

    std::vector<std::pair<size_t, size_t>> noConnections = {};
    std::vector<std::pair<size_t, size_t>> connectionIndices =
        atomPairVectorIndices(atoms, cds::atomPairsConnectedToOtherResidues(pdbResidues));
    writePdbFile(residueIndices, residueTER, writerData, noConnections, outputDir + "glycoprotein_initial");
    resolveOverlapsWithWiggler(assemblyGraphs, assemblyData, glycositeData, proteinResidues_, glycosites_);
    printDihedralAnglesAndOverlapOfGlycosites(assemblyGraphs, assemblyData);
    writerData.atoms.coordinates = cds::atomCoordinates(atoms);
    writePdbFile(residueIndices, residueTER, writerData, noConnections, outputDir + "glycoprotein");
    writerData.atoms.numbers    = cds::serializedNumberVector(atoms.size());
    writerData.residues.numbers = cds::serializedNumberVector(residues.size());
    {
        cds::OffWriterData offData = cds::toOffWriterData(residues);
        offData.atoms.numbers      = writerData.atoms.numbers;
        offData.residues.numbers   = writerData.residues.numbers;
        writeOffFile(offData, outputDir + "glycoprotein");
    }
    writePdbFile(residueIndices, residueTER, writerData, connectionIndices, outputDir + "glycoprotein_serialized");

    for (size_t count = 0; count < settings.number3DStructures; count++)
    {
        resolveOverlapsWithWiggler(assemblyGraphs, assemblyData, glycositeData, proteinResidues_, glycosites_);
        printDihedralAnglesAndOverlapOfGlycosites(assemblyGraphs, assemblyData);
        writerData.atoms.coordinates = cds::atomCoordinates(atoms);
        std::stringstream prefix;
        prefix << count << "_glycoprotein";
        writePdbFile(residueIndices, residueTER, writerData, connectionIndices, outputDir + prefix.str());
    }
}
