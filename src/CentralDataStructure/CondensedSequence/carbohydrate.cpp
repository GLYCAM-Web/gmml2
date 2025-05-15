#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageCreation.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/Assembly/assemblySelection.hpp"
#include "includes/Graph/graphFunctions.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06Functions.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06ResidueNameGenerator.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Graph.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/containerTypes.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"

#include <sstream>
#include <ostream>
#include <cctype>    // isDigit
#include <algorithm> //  std::erase, std::remove

using cdsCondensedSequence::Carbohydrate;

namespace
{
    using cdsCondensedSequence::ParsedResidue;

    std::string getGlycamResidueName(ParsedResidue* residue)
    {
        if (residue->GetType() == cds::ResidueType::Deoxy)
        {
            gmml::log(__LINE__, __FILE__, gmml::WAR,
                      "Bad idea: We asked for Glycam Residue Name of a deoxy type residue (e.g. the 6D of Glc[6D]) "
                      "with name: " +
                          residue->GetResidueName());
            return "";
        }
        std::string linkages = "";
        if (residue->GetType() == cds::ResidueType::Sugar)
        {
            linkages = residue->GetChildLinkagesForGlycamResidueNaming();
        }
        std::string code = GlycamMetadata::Glycam06ResidueNameGenerator(
            linkages, residue->GetPreIsomerModifier(), residue->GetIsomer(), residue->GetResidueName(),
            residue->GetRingType(), residue->GetResidueModifier() + residue->GetRingShape(),
            residue->GetConfiguration());
        return code;
    }

    cds::Coordinate guessCoordinateOfMissingNeighbor(const cds::Atom* centralAtom, double distance)
    {
        if (centralAtom->getNeighbors().size() < 1)
        {
            std::stringstream ss;
            ss << "Error in CreateMissingCoordinateForTetrahedralAtom. centralAtom neighbors is "
               << centralAtom->getNeighbors().size() << " for " << centralAtom->getId();
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
            throw std::runtime_error(ss.str());
        }
        return cds::coordinateOppositeToNeighborAverage(centralAtom->coordinate(),
                                                        cds::atomCoordinates(centralAtom->getNeighbors()), distance);
    }

    // parentAtom (e.g. O of OME), childAtom (e.g. C1 of Gal1-, S1 of SO3)
    void moveConnectedAtomsAccordingToBondLength(cds::Atom* parentAtom, cds::Atom* childAtom)
    {
        double distance      = MolecularMetadata::specificBondLength(parentAtom->getType(), childAtom->getType());
        //  Create an atom c that is will superimpose onto the a atom, bringing b atom with it.
        Coordinate c         = guessCoordinateOfMissingNeighbor(childAtom, distance);
        Coordinate cToParent = parentAtom->coordinate() - c;
        // Figure out which atoms will move
        std::vector<cds::Atom*> atomsToMove;
        atomsToMove.push_back(parentAtom); // add Parent atom so search doesn't go through it.
        cdsSelections::FindConnectedAtoms(atomsToMove, childAtom);
        atomsToMove.erase(atomsToMove.begin()); // delete the parentAtom
        for (auto& atom : atomsToMove)
        {
            atom->setCoordinate(atom->coordinate() + cToParent);
        }
    }

    void derivativeChargeAdjustment(ParsedResidue* parsedResidue)
    {
        std::string adjustAtomName = GlycamMetadata::GetAdjustmentAtom(parsedResidue->getName());
        adjustAtomName             += parsedResidue->GetLinkageName().substr(0, 1);

        cds::Atom* atomToAdjust = parsedResidue->GetParent()->FindAtom(adjustAtomName);
        atomToAdjust->setCharge(atomToAdjust->getCharge() +
                                GlycamMetadata::GetAdjustmentCharge(parsedResidue->getName()));
    }

    cds::Atom* findParentAtom(cds::Residue* parentResidue, cds::Residue* childResidue, const std::string& linkageLabel)
    {
        if (parentResidue->GetType() == cds::ResidueType::Aglycone)
        {
            std::string parentAtomName = GlycamMetadata::GetConnectionAtomForResidue(parentResidue->getName());
            return parentResidue->FindAtom(parentAtomName);
        }
        else if (parentResidue->GetType() == cds::ResidueType::Sugar)
        { // Linkage example: childb1-4parent, it's never parentb1-4child
            size_t linkPosition = 3;
            if (childResidue->GetType() == cds::ResidueType::Derivative)
            { // label will be just a single number.
                linkPosition = 0;
            }
            else if (linkageLabel.size() < 4)
            {
                std::string message =
                    "The deduced linkageLabel is too small:\n" + linkageLabel +
                    ".\nWe require anomer, start atom number, a dash, and connecting atom number. Example:\na1-4";
                gmml::log(__LINE__, __FILE__, gmml::ERR, message);
                throw std::runtime_error(message);
            }
            if (!isdigit(linkageLabel.substr(linkPosition).at(0)))
            {
                std::string message = "Could not convert the last linkage number to an integer: " + linkageLabel;
                gmml::log(__LINE__, __FILE__, gmml::ERR, message);
                throw std::runtime_error(message);
            }
            return cdsSelections::getNonCarbonHeavyAtomNumbered(parentResidue->getAtoms(),
                                                                linkageLabel.substr(linkPosition));
        }
        else
        {
            std::string message = "Error: parent residue: " + parentResidue->getName() +
                                  " isn't either Aglycone or Sugar, and derivatives cannot be parents.";
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
    }

    std::string findChildAtomName(cds::Residue* childResidue, const std::string& linkageLabel)
    {
        std::string childAtomName;
        if (childResidue->GetType() == cds::ResidueType::Derivative)
        {
            std::string glycamNameForResidue =
                getGlycamResidueName(codeUtils::erratic_cast<ParsedResidue*>(childResidue));
            return GlycamMetadata::GetConnectionAtomForResidue(glycamNameForResidue);
        }
        else if (childResidue->GetType() == cds::ResidueType::Sugar)
        {
            std::string childLinkageNumber = linkageLabel.substr(1, 1);
            if (!isdigit(childLinkageNumber.at(0)))
            {
                std::string message = "Could not convert the first linkage number to an integer: " + childLinkageNumber;
                gmml::log(__LINE__, __FILE__, gmml::ERR, message);
                throw std::runtime_error(message);
            }
            return "C" + childLinkageNumber;
        }
        else
        {
            std::string message = "Error: child residue: " + childResidue->getName() +
                                  " is neither derivative or Sugar (aglycones cannot be children)";
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
    }

    void makeDeoxy(cds::Residue* residue, const std::string oxygenNumber)
    { // if oxygenNumber is 6, then C6-O6-H6O becomes C6-Hd
        Atom* hydrogenAtom = residue->FindAtom("H" + oxygenNumber + "O");
        Atom* oxygenAtom   = residue->FindAtom("O" + oxygenNumber);
        Atom* carbonAtom   = residue->FindAtom("C" + oxygenNumber);
        // Add O and H charge to the C atom.
        carbonAtom->setCharge(carbonAtom->getCharge() + oxygenAtom->getCharge() + hydrogenAtom->getCharge());
        // Delete the H of O-H
        residue->deleteAtom(hydrogenAtom);
        // Now transform the Oxygen to a Hd. Easier than deleting O and creating H. Note: this H looks weird in LiteMol
        // as bond length is too long.
        oxygenAtom->setName("Hd");
        oxygenAtom->setType("H1");
        oxygenAtom->setCharge(0.0000);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Completed MakeDeoxy\n");
    }

    void connectAndSetGeometry(cds::Residue* parentResidue, cds::Residue* childResidue)
    {
        std::string linkageLabel = codeUtils::erratic_cast<ParsedResidue*>(childResidue)->GetLinkageName();
        // This is using the new Node<Residue> functionality and the old AtomNode
        // Now go figure out how which Atoms to bond to each other in the residues.
        // Rule: Can't ever have a child aglycone or a parent derivative.
        Atom* parentAtom         = findParentAtom(parentResidue, childResidue, linkageLabel);
        if (parentAtom == nullptr)
        {
            std::string message = "Did not find connection atom in residue: " + residueStringId(parentResidue);
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        // Now get child atom
        std::string childAtomName = findChildAtomName(childResidue, linkageLabel);
        Atom* childAtom           = childResidue->FindAtom(childAtomName);
        if (childAtom == nullptr)
        {
            std::string message = "Did not find child atom named " + childAtomName +
                                  " in child residue: " + residueStringId(childResidue);
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        // Geometry
        moveConnectedAtomsAccordingToBondLength(parentAtom, childAtom);
        //   Now bond the atoms. This could also set distance?, and angle? if passed to function?
        addBond(childAtom, parentAtom); // parentAtom also connected to childAtom. Fancy.
        for (auto& parentAtomNeighbor : parentAtom->getNeighbors())
        {
            if ((parentAtomNeighbor->getName().at(0) != 'H') && (parentAtomNeighbor != childAtom))
            {
                auto matrix =
                    cds::rotationTo(std::array<Coordinate, 3> {parentAtomNeighbor->coordinate(),
                                                               parentAtom->coordinate(), childAtom->coordinate()},
                                    constants::toRadians(constants::DEFAULT_ANGLE));
                std::vector<cds::Atom*> childResidueAtoms = childResidue->mutableAtoms();
                std::vector<cds::Coordinate> coordinates  = cds::atomCoordinates(childResidueAtoms);
                matrix.rotateCoordinates(coordinates);
                cds::setAtomCoordinates(childResidueAtoms, coordinates);
                break;
            }
        }
    }

    std::vector<ParsedResidue*> residuesOrderedByConnectivity(ParsedResidue* terminalResidue)
    {
        std::vector<ParsedResidue*> rawResidues;
        // Go via Graph so order decided by connectivity, depth first traversal:
        glygraph::Graph<cds::Residue> sequenceGraph(terminalResidue);
        for (auto& node : sequenceGraph.getNodes())
        {
            rawResidues.push_back(codeUtils::erratic_cast<ParsedResidue*>(node->getDerivedClass()));
        }
        return rawResidues;
    }

    void setResidueIndices(std::vector<ParsedResidue*> residues)
    {
        unsigned long long linkIndex    = 0; // Convention to start form 0 for linkages.
        unsigned long long residueIndex = 1; // Convention to start from 1 for residues.
        for (auto& residue : residues)
        {
            residue->setIndex(residueIndex);
            residue->setNumber(residueIndex); // ToDo temporary, switch to using number here. Keep index as a gmml
                                              // internal thing, never shown to user.
            ++residueIndex;
            for (auto& edge : residue->getInEdges())
            {
                edge->setIndex(linkIndex);
                ++linkIndex;
            }
        }
        return;
    }

    void createParsedResidues(std::vector<std::unique_ptr<ParsedResidue>>& residuePtrs,
                              const cdsCondensedSequence::SequenceData& sequence)
    {
        size_t residueCount = sequence.residues.name.size();
        residuePtrs.reserve(residueCount);
        for (size_t n = 0; n < residueCount; n++)
        {
            residuePtrs.emplace_back(std::make_unique<ParsedResidue>(cdsCondensedSequence::ParsedResidueComponents {
                sequence.residues.fullString[n],
                sequence.residues.type[n],
                sequence.residues.name[n],
                sequence.residues.linkage[n],
                sequence.residues.ringType[n],
                sequence.residues.configuration[n],
                sequence.residues.isomer[n],
                sequence.residues.preIsomerModifier[n],
                {sequence.residues.ringShape[n], sequence.residues.modifier[n]}
            }));
        }
        for (size_t n = 0; n < sequence.graph.edgeNodes.size(); n++)
        {
            auto& edge = sequence.graph.edgeNodes[n];
            residuePtrs[edge[1]].get()->addParent(sequence.edges.names[n], residuePtrs[edge[0]].get());
        }
    }

    void sortResidueEdges(std::vector<std::unique_ptr<ParsedResidue>>& residuePtrs)
    {
        for (auto& residue : residuePtrs)
        {
            residue.get()->sortOutEdgesBySourceTObjectComparator();
        }
    }

    void initialWiggleLinkage(const codeUtils::SparseVector<double>& elementRadii,
                              const GlycamMetadata::DihedralAngleDataTable& metadataTable, cds::Molecule* molecule,
                              cds::Residue* residue, cds::ResidueLinkage& linkage,
                              const cds::AngleSearchSettings& searchSettings)
    {
        // GREEDY: taken care of, but note that the atoms that move in RotatableDihedral class need to be updated after
        // more residues are added.
        auto shapePreference = cds::firstRotamerOnly(
            linkage, cds::defaultShapePreference(metadataTable, linkage.rotamerType, linkage.dihedralMetadata));
        cds::setShapeToPreference(linkage, shapePreference);
        auto searchPreference               = cds::angleSearchPreference(searchSettings.deviation, shapePreference);
        const cds::GraphIndexData graphData = cds::toIndexData({molecule});
        const assembly::Graph graph         = cds::createCompleteAssemblyGraph(graphData);
        size_t residueIndex                 = codeUtils::indexOf(graphData.objects.residues, residue);
        std::vector<bool> reachable =
            graph::reachableNodes(graph.residues, std::vector<bool>(graph.indices.residueCount, false), residueIndex);
        const assembly::Selection selection = assembly::selectByResidues(graph, reachable);
        const assembly::Bounds bounds       = toAssemblyBounds(elementRadii, graphData, graph);
        std::vector<std::array<std::vector<bool>, 2>> residueAtomsCloseToEdge =
            assembly::atomsCloseToResidueEdges(graph);
        assembly::Bounds newBounds = cds::simpleWiggleCurrentRotamers(
            GlycamMetadata::dihedralAngleDataTable(), MolecularMetadata::potentialTable(), constants::overlapTolerance,
            searchSettings.angles, linkage.rotatableDihedrals, linkage.dihedralMetadata, searchPreference,
            graphData.objects, graph, selection, bounds, residueAtomsCloseToEdge);
        for (size_t n = 0; n < graph.indices.atomCount; n++)
        {
            graphData.objects.atoms[n]->setCoordinate(newBounds.atoms[n].center);
        }
    }

    // Gonna choke on cyclic glycans. Add a check for IsVisited when that is required.
    void depthFirstSetConnectivityAndGeometry(cds::Molecule* molecule,
                                              std::vector<cds::ResidueLinkage>& glycosidicLinkages,
                                              const cds::AngleSearchSettings& searchSettings,
                                              const codeUtils::SparseVector<double>& elementRadii,
                                              const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                              cds::Residue* currentParent)
    {
        // MolecularModeling::ResidueVector neighbors = to_this_residue2->GetNode()->GetResidueNeighbors();

        // Additional code to sort neighbors by lowest index.
        // Only required so that numbers match those assigned in condensed sequence class
        // Should not be done this way, need a generic graph structure and then to centralize everything.
        // MolecularModeling::ResidueVector neighbors =
        // selection::SortResidueNeighborsByAcendingConnectionAtomNumber(to_this_residue2->GetNode()->GetResidueNodeConnectingAtoms());
        // End addtional sorting code.
        // Breath first code
        // for(auto &neighbor : neighbors)
        // {
        //     if(neighbor->GetIndex() != from_this_residue1->GetIndex()) // If not the previous residue
        //     {
        //         residue_linkages->emplace_back(neighbor, to_this_residue2);
        //     }
        // }
        // End Breath first code
        for (auto& child : currentParent->getChildren())
        {
            connectAndSetGeometry(currentParent, child);
            cds::ResidueLink link = cds::findResidueLink({child, currentParent});
            cds::ResidueLinkage& linkage =
                glycosidicLinkages.emplace_back(cds::createResidueLinkage(metadataTable, link));
            initialWiggleLinkage(elementRadii, metadataTable, molecule, child, linkage, searchSettings);
            depthFirstSetConnectivityAndGeometry(molecule, glycosidicLinkages, searchSettings, elementRadii,
                                                 metadataTable, child);
        }
    }
} // namespace

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Carbohydrate::Carbohydrate(const cdsParameters::ParameterManager& parameterManager,
                           const codeUtils::SparseVector<double>& elementRadii, const std::string& inputSequence)
    : cds::Molecule()
{
    {
        std::vector<std::unique_ptr<ParsedResidue>> residuePtrs;
        cdsCondensedSequence::SequenceData sequenceData = reordered(parseSequence(inputSequence));
        size_t residueCount                             = sequenceData.residues.name.size();
        createParsedResidues(residuePtrs, sequenceData);
        sortResidueEdges(residuePtrs);

        for (auto& residue : codeUtils::pointerToUniqueVector(residuePtrs))
        { // Move atoms from prep file into parsedResidues.
            if (residue->GetType() != cds::ResidueType::Deoxy)
            {
                cdsParameters::createAtomsForResidue(parameterManager, residue, getGlycamResidueName(residue));
                if (residue->GetType() == cds::ResidueType::Derivative)
                { // Deal with adjusting charges for derivatives
                    derivativeChargeAdjustment(residue);
                }
            }
        }
        for (size_t n = residueCount - 1; n < residueCount; n--)
        { // Apply any deoxy
            if (sequenceData.residues.type[n] == cds::ResidueType::Deoxy)
            {
                ParsedResidue* deoxyResidue         = residuePtrs[n].get();
                ParsedResidue* residueToBeDeoxified = deoxyResidue->GetParent();
                makeDeoxy(residueToBeDeoxified, deoxyResidue->GetLink());
                residuePtrs.erase(residuePtrs.begin() + n);
            }
        }
        setResidueIndices(
            residuesOrderedByConnectivity(residuePtrs[0].get())); // For reporting residue index numbers to the user
        for (auto& res : residuePtrs)
        {
            this->addResidue(std::move(res));
        }
    }
    this->setName("CONDENSEDSEQUENCE");
    // Have atom numbers go from 1 to number of atoms. Note this should be after deleting atoms due to deoxy

    cds::serializeNumbers(this->getAtoms());
    auto searchSettings                                         = defaultSearchSettings;
    // Set 3D structure
    const GlycamMetadata::DihedralAngleDataTable& metadataTable = GlycamMetadata::dihedralAngleDataTable();
    depthFirstSetConnectivityAndGeometry(this, glycosidicLinkages_, searchSettings, elementRadii, metadataTable,
                                         getResidues().front()); // recurve start with terminal
    // Re-numbering is a hack as indices have global scope and two instances give too high numbers.
    unsigned int linkageIndex = 0;
    // Linkages should be Edges to avoid this as they already get renumbered above.
    for (auto& linkage : glycosidicLinkages_) // These will exist on the vector in order of edge connectivity set above.
    { // Greedy first means the atoms-to-move needs to be updated for every linkage:
        linkage.index = linkageIndex++;
        cds::determineAtomsThatMove(linkage.rotatableDihedrals);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Final carbohydrate overlap resolution starting.");
    this->ResolveOverlaps(elementRadii, metadataTable, searchSettings);
    gmml::log(__LINE__, __FILE__, gmml::INF,
              "Final carbohydrate overlap resolution finished. Returning from carbohydrate ctor");
    return;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void Carbohydrate::deleteResidue(cds::Residue* byeBye)
{
    if (!glycosidicLinkages_.empty())
    {
        throw std::runtime_error("linkages should not have been created at this point");
    }
    cds::Molecule::deleteResidue(byeBye);
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
// std::string fileOutputDirectory = "unspecified", std::string fileType = "PDB", std::string outputFileNaming =
// "structure"
void Carbohydrate::Generate3DStructureFiles(const std::string& fileOutputDirectory, const std::string& outputFileNaming,
                                            const std::vector<std::string>& headerLines)
{ // ToDo exception handling in centralized function for writing pdb/off
    // Build the filename and path, add appropriate suffix later
    try
    {
        std::string PathAndFileName;
        if (fileOutputDirectory == "unspecified") // "unspecified" is the default
        {
            PathAndFileName += "./" + outputFileNaming;
        }
        else
        {
            PathAndFileName += fileOutputDirectory + "/" + outputFileNaming;
        }
        cds::GraphIndexData indices = cds::toIndexData({this});
        codeUtils::writeToFile(PathAndFileName + ".pdb",
                               [&](std::ostream& stream)
                               {
                                   WritePdb(stream, indices, headerLines);
                               });
        codeUtils::writeToFile(PathAndFileName + ".off",
                               [&](std::ostream& stream)
                               {
                                   WriteOff(stream, getName(), indices);
                               });
    }
    catch (const std::string& exceptionMessage)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "carbohydrate class caught this exception message: " + exceptionMessage);
        throw exceptionMessage;
    }
    catch (const std::runtime_error& error)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, error.what());
        throw error;
    }
    catch (...)
    {
        std::string message = "carbohydrate class caught a throw type that was not anticipated. Pretty please report "
                              "how you got to this to glycam@gmail.com.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
}

std::string Carbohydrate::GetNumberOfShapes(bool likelyShapesOnly) const
{
    if (this->CountShapes(likelyShapesOnly) > 4294967296)
    {
        return ">2^32";
    }
    return std::to_string(this->CountShapes(likelyShapesOnly));
}

unsigned long int Carbohydrate::CountShapes(bool likelyShapesOnly) const
{
    const GlycamMetadata::DihedralAngleDataTable& metadataTable = GlycamMetadata::dihedralAngleDataTable();
    unsigned long long int numberOfShapes                       = 1;
    for (auto& linkage : cds::nonDerivativeResidueLinkages(glycosidicLinkages_))
    {
        auto rotamerType = linkage.rotamerType;
        auto& metadata   = linkage.dihedralMetadata;
        numberOfShapes   *= likelyShapesOnly ? cds::numberOfLikelyShapes(metadataTable, rotamerType, metadata)
                                             : cds::numberOfShapes(metadataTable, rotamerType, metadata);
    }
    return numberOfShapes;
}

//////////////////////////////////////////////////////////
//                  PRIVATE FUNCTIONS                   //
//////////////////////////////////////////////////////////

void Carbohydrate::ResolveOverlaps(const codeUtils::SparseVector<double>& elementRadii,
                                   const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                   const cds::AngleSearchSettings& searchSettings)
{
    const cds::GraphIndexData graphData = cds::toIndexData({this});
    const assembly::Graph graph         = cds::createCompleteAssemblyGraph(graphData);
    const assembly::Selection selection = assembly::selectAll(graph);
    assembly::Bounds bounds             = cds::toAssemblyBounds(elementRadii, graphData, graph);
    std::vector<std::array<std::vector<bool>, 2>> residueAtomsCloseToEdge = assembly::atomsCloseToResidueEdges(graph);
    // wiggle twice for nicer structures
    for (size_t n = 0; n < 2; n++)
    {
        for (auto& linkage : glycosidicLinkages_)
        {
            auto preference = cds::angleSearchPreference(
                searchSettings.deviation,
                cds::currentRotamerOnly(linkage, cds::defaultShapePreference(metadataTable, linkage.rotamerType,
                                                                             linkage.dihedralMetadata)));
            bounds = cds::simpleWiggleCurrentRotamers(
                metadataTable, MolecularMetadata::potentialTable(), constants::overlapTolerance, searchSettings.angles,
                linkage.rotatableDihedrals, linkage.dihedralMetadata, preference, graphData.objects, graph, selection,
                bounds, residueAtomsCloseToEdge);
        }
    }
    for (size_t n = 0; n < graph.indices.atomCount; n++)
    {
        graphData.objects.atoms[n]->setCoordinate(bounds.atoms[n].center);
    }
    return;
}
