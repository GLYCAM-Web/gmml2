#include "include/carbohydrate/carbohydrate.hpp"

#include "include/CentralDataStructure/Selections/atomSelections.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/CentralDataStructure/molecule.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/CentralDataStructure/residueLinkage/residueLinkageCreation.hpp"
#include "include/CentralDataStructure/residueLinkage/residueLinkageFunctions.hpp"
#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblySelection.hpp"
#include "include/carbohydrate/dihedralAngleSearch.hpp"
#include "include/carbohydrate/dihedralShape.hpp"
#include "include/carbohydrate/offWriter.hpp"
#include "include/carbohydrate/pdbWriter.hpp"
#include "include/fileType/off/offFileWriter.hpp"
#include "include/fileType/pdb/pdbFileWriter.hpp"
#include "include/fileType/pdb/pdbResidue.hpp"
#include "include/geometry/matrix.hpp"
#include "include/geometry/measurements.hpp"
#include "include/geometry/orientation.hpp"
#include "include/graph/graphFunctions.hpp"
#include "include/metadata/atomicBonds.hpp"
#include "include/metadata/glycam06Functions.hpp"
#include "include/metadata/glycam06ResidueNameGenerator.hpp"
#include "include/preprocess/parameterManager.hpp"
#include "include/sequence/parsedResidue.hpp"
#include "include/sequence/sequenceTypes.hpp"
#include "include/sequence/sequenceUtil.hpp"
#include "include/templateGraph/Graph.hpp"
#include "include/util/casting.hpp"
#include "include/util/constants.hpp"
#include "include/util/containerTypes.hpp"
#include "include/util/containers.hpp"
#include "include/util/files.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <algorithm> //  std::erase, std::remove
#include <cctype>    // isDigit
#include <memory>
#include <ostream>
#include <sstream>

namespace gmml
{
    namespace
    {
        template<class T> std::vector<T*> pointerToUniqueVector(std::vector<std::unique_ptr<T>>& vec)
        {
            std::vector<T*> result;
            result.reserve(vec.size());
            for (auto& a : vec)
            {
                result.push_back(a.get());
            }
            return result;
        }

        std::string childLinkagesForGlycamResidueNaming(const ParsedResidue* residue)
        {
            const std::vector<ParsedResidue*> children = residue->GetChildren();
            std::vector<std::string> nonDeoxyLinkNames;
            nonDeoxyLinkNames.reserve(children.size());
            for (auto child : children)
            {
                if (child->GetType() != ResidueType::Deoxy)
                { // For glycam residue name, e.g. 2YB, do not want deoxy linkages to impact the residue name.
                    nonDeoxyLinkNames.push_back(child->GetLink());
                }
            }
            if (nonDeoxyLinkNames.empty())
            {
                return "Terminal";
            }
            else
            {
                return util::join(",", nonDeoxyLinkNames);
            }
        }

        std::string getGlycamResidueName(ParsedResidue* residue)
        {
            if (residue->GetType() == ResidueType::Deoxy)
            {
                util::log(
                    __LINE__,
                    __FILE__,
                    util::WAR,
                    "Bad idea: We asked for Glycam Residue Name of a deoxy type residue (e.g. the 6D of Glc[6D]) "
                    "with name: " +
                        residue->GetResidueName());
                return "";
            }
            std::string linkages = "";
            if (residue->GetType() == ResidueType::Sugar)
            {
                linkages = childLinkagesForGlycamResidueNaming(residue);
            }
            std::string code = metadata::Glycam06ResidueNameGenerator(
                linkages,
                residue->GetPreIsomerModifier(),
                residue->GetIsomer(),
                residue->GetResidueName(),
                residue->GetRingType(),
                residue->GetResidueModifier() + residue->GetRingShape(),
                residue->GetConfiguration());
            return code;
        }

        Coordinate guessCoordinateOfMissingNeighbor(const Atom* centralAtom, double distance)
        {
            if (centralAtom->getNeighbors().size() < 1)
            {
                std::stringstream ss;
                ss << "Error in CreateMissingCoordinateForTetrahedralAtom. centralAtom neighbors is "
                   << centralAtom->getNeighbors().size() << " for " << centralAtom->getId();
                util::log(__LINE__, __FILE__, util::ERR, ss.str());
                throw std::runtime_error(ss.str());
            }
            return coordinateOppositeToNeighborAverage(
                centralAtom->coordinate(), atomCoordinates(centralAtom->getNeighbors()), distance);
        }

        // parentAtom (e.g. O of OME), childAtom (e.g. C1 of Gal1-, S1 of SO3)
        void moveConnectedAtomsAccordingToBondLength(Atom* parentAtom, Atom* childAtom)
        {
            double distance = specificBondLength(parentAtom->getType(), childAtom->getType());
            //  Create an atom c that is will superimpose onto the a atom, bringing b atom with it.
            Coordinate c = guessCoordinateOfMissingNeighbor(childAtom, distance);
            Coordinate cToParent = parentAtom->coordinate() - c;
            // Figure out which atoms will move
            std::vector<Atom*> atomsToMove;
            atomsToMove.push_back(parentAtom); // add Parent atom so search doesn't go through it.
            FindConnectedAtoms(atomsToMove, childAtom);
            atomsToMove.erase(atomsToMove.begin()); // delete the parentAtom
            for (auto& atom : atomsToMove)
            {
                atom->setCoordinate(atom->coordinate() + cToParent);
            }
        }

        void derivativeChargeAdjustment(ParsedResidue* parsedResidue)
        {
            std::string adjustAtomName = metadata::GetAdjustmentAtom(parsedResidue->getName());
            adjustAtomName += parsedResidue->GetLinkageName().substr(0, 1);

            Atom* atomToAdjust = parsedResidue->GetParent()->FindAtom(adjustAtomName);
            atomToAdjust->setCharge(
                atomToAdjust->getCharge() + metadata::GetAdjustmentCharge(parsedResidue->getName()));
        }

        Atom* findParentAtom(Residue* parentResidue, Residue* childResidue, const std::string& linkageLabel)
        {
            if (parentResidue->GetType() == ResidueType::Aglycone)
            {
                std::string parentAtomName = metadata::GetConnectionAtomForResidue(parentResidue->getName());
                return parentResidue->FindAtom(parentAtomName);
            }
            else if (parentResidue->GetType() == ResidueType::Sugar)
            { // Linkage example: childb1-4parent, it's never parentb1-4child
                size_t linkPosition = 3;
                if (childResidue->GetType() == ResidueType::Derivative)
                { // label will be just a single number.
                    linkPosition = 0;
                }
                else if (linkageLabel.size() < 4)
                {
                    std::string message =
                        "The deduced linkageLabel is too small:\n" + linkageLabel +
                        ".\nWe require anomer, start atom number, a dash, and connecting atom number. Example:\na1-4";
                    util::log(__LINE__, __FILE__, util::ERR, message);
                    throw std::runtime_error(message);
                }
                if (!isdigit(linkageLabel.substr(linkPosition).at(0)))
                {
                    std::string message = "Could not convert the last linkage number to an integer: " + linkageLabel;
                    util::log(__LINE__, __FILE__, util::ERR, message);
                    throw std::runtime_error(message);
                }
                return getNonCarbonHeavyAtomNumbered(parentResidue->getAtoms(), linkageLabel.substr(linkPosition));
            }
            else
            {
                std::string message = "Error: parent residue: " + parentResidue->getName() +
                                      " isn't either Aglycone or Sugar, and derivatives cannot be parents.";
                util::log(__LINE__, __FILE__, util::ERR, message);
                throw std::runtime_error(message);
            }
        }

        std::string findChildAtomName(Residue* childResidue, const std::string& linkageLabel)
        {
            std::string childAtomName;
            if (childResidue->GetType() == ResidueType::Derivative)
            {
                std::string glycamNameForResidue =
                    getGlycamResidueName(util::erratic_cast<ParsedResidue*>(childResidue));
                return metadata::GetConnectionAtomForResidue(glycamNameForResidue);
            }
            else if (childResidue->GetType() == ResidueType::Sugar)
            {
                std::string childLinkageNumber = linkageLabel.substr(1, 1);
                if (!isdigit(childLinkageNumber.at(0)))
                {
                    std::string message =
                        "Could not convert the first linkage number to an integer: " + childLinkageNumber;
                    util::log(__LINE__, __FILE__, util::ERR, message);
                    throw std::runtime_error(message);
                }
                return "C" + childLinkageNumber;
            }
            else
            {
                std::string message = "Error: child residue: " + childResidue->getName() +
                                      " is neither derivative or Sugar (aglycones cannot be children)";
                util::log(__LINE__, __FILE__, util::ERR, message);
                throw std::runtime_error(message);
            }
        }

        void makeDeoxy(Residue* residue, const std::string oxygenNumber)
        { // if oxygenNumber is 6, then C6-O6-H6O becomes C6-Hd
            Atom* hydrogenAtom = residue->FindAtom("H" + oxygenNumber + "O");
            Atom* oxygenAtom = residue->FindAtom("O" + oxygenNumber);
            Atom* carbonAtom = residue->FindAtom("C" + oxygenNumber);
            // Add O and H charge to the C atom.
            carbonAtom->setCharge(carbonAtom->getCharge() + oxygenAtom->getCharge() + hydrogenAtom->getCharge());
            // Delete the H of O-H
            residue->deleteAtom(hydrogenAtom);
            // Now transform the Oxygen to a Hd. Easier than deleting O and creating H
            Coordinate carbonPos = carbonAtom->coordinate();
            Coordinate direction = normal(oxygenAtom->coordinate() - carbonPos);
            oxygenAtom->setCoordinate(carbonPos + scaleBy(deoxyHydrogenBondLength(), direction));
            oxygenAtom->setName("Hd");
            oxygenAtom->setType("H1");
            oxygenAtom->setCharge(0.0000);
            util::log(__LINE__, __FILE__, util::INF, "Completed MakeDeoxy\n");
        }

        void connectAndSetGeometry(Residue* parentResidue, Residue* childResidue)
        {
            std::string linkageLabel = util::erratic_cast<ParsedResidue*>(childResidue)->GetLinkageName();
            // This is using the new Node<Residue> functionality and the old AtomNode
            // Now go figure out how which Atoms to bond to each other in the residues.
            // Rule: Can't ever have a child aglycone or a parent derivative.
            Atom* parentAtom = findParentAtom(parentResidue, childResidue, linkageLabel);
            if (parentAtom == nullptr)
            {
                std::string message = "Did not find connection atom in residue: " + residueStringId(parentResidue);
                util::log(__LINE__, __FILE__, util::ERR, message);
                throw std::runtime_error(message);
            }
            // Now get child atom
            std::string childAtomName = findChildAtomName(childResidue, linkageLabel);
            Atom* childAtom = childResidue->FindAtom(childAtomName);
            if (childAtom == nullptr)
            {
                std::string message = "Did not find child atom named " + childAtomName +
                                      " in child residue: " + residueStringId(childResidue);
                util::log(__LINE__, __FILE__, util::ERR, message);
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
                    auto matrix = rotationTo(
                        std::array<Coordinate, 3> {
                            parentAtomNeighbor->coordinate(), parentAtom->coordinate(), childAtom->coordinate()},
                        constants::toRadians(constants::DEFAULT_ANGLE));
                    std::vector<Atom*> childResidueAtoms = childResidue->mutableAtoms();
                    std::vector<Coordinate> coordinates = atomCoordinates(childResidueAtoms);
                    setAtomCoordinates(childResidueAtoms, transform(matrix, coordinates));
                    break;
                }
            }
        }

        std::vector<ParsedResidue*> residuesOrderedByConnectivity(ParsedResidue* terminalResidue)
        {
            std::vector<ParsedResidue*> rawResidues;
            // Go via Graph so order decided by connectivity, depth first traversal:
            glygraph::Graph<Residue> sequenceGraph(terminalResidue);
            for (auto& node : sequenceGraph.getNodes())
            {
                rawResidues.push_back(util::erratic_cast<ParsedResidue*>(node->getDerivedClass()));
            }
            return rawResidues;
        }

        void setResidueIndices(std::vector<ParsedResidue*> residues)
        {
            unsigned long long linkIndex = 0;    // Convention to start form 0 for linkages.
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

        void createParsedResidues(
            std::vector<std::unique_ptr<ParsedResidue>>& residuePtrs,
            std::vector<size_t>& indices,
            const sequence::SequenceData& sequence)
        {
            size_t residueCount = nodeCount(sequence.graph);
            residuePtrs.reserve(residueCount);
            indices.reserve(residueCount);
            std::vector<size_t> newIndices;
            newIndices.reserve(residueCount);
            size_t index = 0;
            for (size_t n = 0; n < residueCount; n++)
            {
                if (sequence.graph.nodes.alive[n])
                {
                    size_t edge = parentEdge(sequence, n);
                    bool hasParent = edge < edgeCount(sequence.graph);
                    residuePtrs.emplace_back(std::make_unique<ParsedResidue>(sequence::ParsedResidueComponents {
                        sequence.residues.fullString[n],
                        sequence.residues.type[n],
                        sequence.residues.name[n],
                        (hasParent ? edgeLinkage(sequence, edge) : ""),
                        "",
                        sequence.residues.ringType[n],
                        sequence.residues.configuration[n],
                        sequence.residues.isomer[n],
                        sequence.residues.preIsomerModifier[n],
                        sequence.residues.ringShape[n],
                        sequence.residues.modifier[n]}));
                    indices.push_back(n);
                    newIndices.push_back(index);
                    index++;
                }
                else
                {
                    newIndices.push_back(residueCount);
                }
            }
            for (size_t n = 0; n < edgeCount(sequence.graph); n++)
            {
                auto& edge = sequence.graph.edges.nodes[n];
                if (sequence.graph.nodes.alive[edge[1]] && !sequence.graph.edges.alive[edge[0]])
                {
                    throw std::runtime_error(
                        "Error: required parent residue removed by random chance. Check sequence definition");
                }
                if (graph::edgeAlive(sequence.graph, n))
                {
                    size_t childId = newIndices[edge[1]];
                    size_t parentId = newIndices[edge[0]];
                    residuePtrs[parentId].get()->addChild(sequence.edges.names[n], residuePtrs[childId].get());
                }
            }
        }

        std::vector<DihedralIndices> linkageDihedralIndices(
            const GraphIndexData& graphData, const ResidueLinkage& linkage)
        {
            std::vector<DihedralIndices> result;
            result.reserve(linkage.rotatableBonds.size());
            for (auto& bond : linkage.rotatableBonds)
            {
                auto index = [&](size_t n) { return util::indexOf(graphData.objects.atoms, bond.dihedralAtoms[n]); };
                result.push_back({
                    {index(0), index(1), index(2), index(3)},
                    util::indicesOf(graphData.objects.atoms, bond.movingAtoms),
                    bond.currentMetadataIndex
                });
            }
            return result;
        }

        void initialWiggleLinkage(
            const util::SparseVector<double>& elementRadii,
            const DihedralAngleDataTable& dihedralAngleData,
            Molecule& molecule,
            Residue* residue,
            ResidueLinkage& linkage,
            const AngleSearchSettings& searchSettings)
        {
            // GREEDY: taken care of, but note that the atoms that move in RotatableDihedral class need to be updated
            // after more residues are added.
            auto shapePreference = firstRotamerOnly(
                linkage, defaultShapePreference(dihedralAngleData, linkage.rotamerType, linkage.dihedralMetadata));
            setShapeToPreference(linkage, shapePreference);
            auto searchPreference = angleSearchPreference(searchSettings.deviation, shapePreference);
            const GraphIndexData graphData = toIndexData({&molecule});
            const assembly::Graph graph = createCompleteAssemblyGraph(graphData);
            size_t residueIndex = util::indexOf(graphData.objects.residues, residue);
            std::vector<bool> reachable = graph::reachableNodes(
                graph.residues, std::vector<bool>(residueCount(graph.source), false), residueIndex);
            const assembly::Selection selection = assembly::selectByResidues(graph, reachable);
            const std::vector<Sphere> atomBounds = assembly::toAtomBounds(
                elementRadii, atomElements(graphData.objects.atoms), atomCoordinates(graphData.objects.atoms));
            const assembly::Bounds bounds = toAssemblyBounds(graph, atomBounds);
            std::vector<std::array<std::vector<bool>, 2>> residueAtomsCloseToEdge =
                assembly::atomsCloseToResidueEdges(graph);
            std::vector<bool> foundElements = gmml::foundElements(atomElements(graphData.objects.atoms));
            const PotentialTable potentials = potentialTable(elementRadii, foundElements);
            assembly::Bounds newBounds = simpleWiggleCurrentRotamers(
                dihedralAngleData,
                potentials,
                searchSettings.angles,
                linkageDihedralIndices(graphData, linkage),
                linkage.dihedralMetadata,
                searchPreference,
                atomElements(graphData.objects.atoms),
                graph,
                selection,
                bounds,
                residueAtomsCloseToEdge);
            for (size_t n = 0; n < atomCount(graph.source); n++)
            {
                graphData.objects.atoms[n]->setCoordinate(newBounds.atoms[n].center);
            }
        }

        // Gonna choke on cyclic glycans. Add a check for IsVisited when that is required.
        void depthFirstSetConnectivityAndGeometry(
            Molecule& molecule,
            std::vector<ResidueLinkage>& glycosidicLinkages,
            const AngleSearchSettings& searchSettings,
            const util::SparseVector<double>& elementRadii,
            const DihedralAngleDataTable& dihedralAngleData,
            Residue* currentParent)
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
                ResidueLink link = findResidueLink({child, currentParent});
                ResidueLinkage& linkage =
                    glycosidicLinkages.emplace_back(createResidueLinkage(dihedralAngleData, link));
                initialWiggleLinkage(elementRadii, dihedralAngleData, molecule, child, linkage, searchSettings);
                depthFirstSetConnectivityAndGeometry(
                    molecule, glycosidicLinkages, searchSettings, elementRadii, dihedralAngleData, child);
            }
        }
    } // namespace

    void initializeCarbohydrate(
        Molecule& molecule,
        std::vector<ResidueLinkage>& glycosidicLinkages,
        const preprocess::ParameterManager& parameterManager,
        const util::SparseVector<double>& elementRadii,
        const DihedralAngleDataTable& dihedralAngleData,
        const sequence::SequenceData& sequence)
    {
        {
            std::vector<std::unique_ptr<ParsedResidue>> residuePtrs;
            std::vector<size_t> indices;
            createParsedResidues(residuePtrs, indices, sequence);
            size_t residueCount = residuePtrs.size();

            for (auto& residue : pointerToUniqueVector(residuePtrs))
            { // Move atoms from prep file into parsedResidues.
                if (residue->GetType() != ResidueType::Deoxy)
                {
                    createAtomsForResidue(parameterManager, residue, getGlycamResidueName(residue));
                    if (residue->GetType() == ResidueType::Derivative)
                    { // Deal with adjusting charges for derivatives
                        derivativeChargeAdjustment(residue);
                    }
                }
            }
            const sequence::ResidueData& rD = sequence.residues;
            std::vector<ResidueType> residueTypes = util::indicesToValues(rD.type, indices);
            for (size_t n = residueCount - 1; n < residueCount; n--)
            { // Apply any deoxy and set residue attributes
                size_t k = indices[n];
                size_t edge = parentEdge(sequence, k);
                std::string linkage = (edge < edgeCount(sequence.graph)) ? edgeLinkage(sequence, edge) : "";
                if (residueTypes[n] == ResidueType::Deoxy)
                {
                    ParsedResidue* deoxyResidue = residuePtrs[n].get();
                    ParsedResidue* residueToBeDeoxified = deoxyResidue->GetParent();
                    makeDeoxy(residueToBeDeoxified, linkage);
                    residuePtrs.erase(residuePtrs.begin() + n);
                }
                else
                {
                    std::string glycamCode = getGlycamResidueName(residuePtrs[n].get());
                    ResidueAttributes ra {
                        rD.type[k],
                        rD.name[k],
                        glycamCode,
                        linkage,
                        rD.ringType[k],
                        rD.configuration[k],
                        rD.isomer[k],
                        rD.preIsomerModifier[k],
                        rD.modifier[k],
                        rD.isInternal[k]};
                    residuePtrs[n]->setAttributes(ra);
                }
            }
            setResidueIndices(
                residuesOrderedByConnectivity(residuePtrs[0].get())); // For reporting residue index numbers to the user
            for (auto& res : residuePtrs)
            {
                molecule.addResidue(std::move(res));
            }
        }
        molecule.setName("CONDENSEDSEQUENCE");
        // Have atom numbers go from 1 to number of atoms. Note this should be after deleting atoms due to deoxy

        serializeNumbers(molecule.getAtoms());
        auto searchSettings = defaultSearchSettings;
        // Set 3D structure
        depthFirstSetConnectivityAndGeometry(
            molecule,
            glycosidicLinkages,
            searchSettings,
            elementRadii,
            dihedralAngleData,
            molecule.getResidues().front()); // recurve start with terminal
        // Re-numbering is a hack as indices have global scope and two instances give too high numbers.
        unsigned int linkageIndex = 0;
        // Linkages should be Edges to avoid this as they already get renumbered above.
        for (auto& linkage :
             glycosidicLinkages) // These will exist on the vector in order of edge connectivity set above.
        {                        // Greedy first means the atoms-to-move needs to be updated for every linkage:
            linkage.index = linkageIndex++;
            determineAtomsThatMove(linkage.rotatableBonds);
        }
        util::log(__LINE__, __FILE__, util::INF, "Final carbohydrate overlap resolution starting.");
        resolveOverlaps(molecule, glycosidicLinkages, elementRadii, dihedralAngleData, searchSettings);
        util::log(
            __LINE__,
            __FILE__,
            util::INF,
            "Final carbohydrate overlap resolution finished. Returning from carbohydrate initializer");
    }

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    // std::string fileOutputDirectory = "unspecified", std::string fileType = "PDB", std::string outputFileNaming =
    // "structure"
    void generate3DStructureFiles(
        const carbohydrate::CarbohydrateData& data,
        const std::string& name,
        const std::string& fileOutputDirectory,
        const std::string& outputFileNaming,
        const std::vector<std::string>& headerLines)
    { // ToDo exception handling in centralized function for writing pdb/off
        // Build the filename and path, add appropriate suffix later
        assembly::Graph graph = createVisibleAssemblyGraph(data);
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
            util::writeToFile(
                PathAndFileName + ".pdb", [&](std::ostream& stream) { writePdb(stream, data, graph, headerLines); });
            util::writeToFile(
                PathAndFileName + ".off",
                [&](std::ostream& stream)
                {
                    off::OffFileData offData = toOffFileData(data, graph);
                    off::writeResiduesTogether(stream, offData, graph, indices(graph.residues.nodes), name);
                });
        }
        catch (const std::string& exceptionMessage)
        {
            util::log(
                __LINE__, __FILE__, util::ERR, "carbohydrate class caught this exception message: " + exceptionMessage);
            throw exceptionMessage;
        }
        catch (const std::runtime_error& error)
        {
            util::log(__LINE__, __FILE__, util::ERR, error.what());
            throw error;
        }
        catch (...)
        {
            std::string message =
                "carbohydrate class caught a throw type that was not anticipated. Pretty please report "
                "how you got to this to glycam@gmail.com.";
            util::log(__LINE__, __FILE__, util::ERR, message);
            throw std::runtime_error(message);
        }
    }

    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////

    void resolveOverlaps(
        Molecule& molecule,
        std::vector<ResidueLinkage>& glycosidicLinkages,
        const util::SparseVector<double>& elementRadii,
        const DihedralAngleDataTable& dihedralAngleData,
        const AngleSearchSettings& searchSettings)
    {
        const GraphIndexData graphData = toIndexData({&molecule});
        const assembly::Graph graph = createCompleteAssemblyGraph(graphData);
        const assembly::Selection selection = selectAll(graph);
        const std::vector<Sphere> atomBounds = assembly::toAtomBounds(
            elementRadii, atomElements(graphData.objects.atoms), atomCoordinates(graphData.objects.atoms));
        assembly::Bounds bounds = toAssemblyBounds(graph, atomBounds);
        std::vector<Element> atomElements = gmml::atomElements(graphData.objects.atoms);
        std::vector<bool> foundElements = gmml::foundElements(atomElements);
        const PotentialTable potentials = potentialTable(elementRadii, foundElements);
        std::vector<std::array<std::vector<bool>, 2>> residueAtomsCloseToEdge =
            assembly::atomsCloseToResidueEdges(graph);
        // wiggle twice for nicer structures
        for (size_t n = 0; n < 2; n++)
        {
            for (auto& linkage : glycosidicLinkages)
            {
                auto preference = angleSearchPreference(
                    searchSettings.deviation,
                    currentRotamerOnly(
                        linkage,
                        defaultShapePreference(dihedralAngleData, linkage.rotamerType, linkage.dihedralMetadata)));
                bounds = simpleWiggleCurrentRotamers(
                    dihedralAngleData,
                    potentials,
                    searchSettings.angles,
                    linkageDihedralIndices(graphData, linkage),
                    linkage.dihedralMetadata,
                    preference,
                    atomElements,
                    graph,
                    selection,
                    bounds,
                    residueAtomsCloseToEdge);
            }
        }
        for (size_t n = 0; n < atomCount(graph.source); n++)
        {
            graphData.objects.atoms[n]->setCoordinate(bounds.atoms[n].center);
        }
        return;
    }

    carbohydrate::CarbohydrateData structured(const std::vector<Molecule*>& molecules)
    {
        GraphIndexData indices = toIndexData(molecules);
        assembly::Graph graph = createCompleteAssemblyGraph(indices);
        graph::Graph atomGraph = graph.atoms;
        std::vector<Atom*>& atoms = indices.objects.atoms;
        std::vector<Residue*>& residues = indices.objects.residues;
        carbohydrate::AtomData atomData {
            atomNames(atoms),
            atomTypes(atoms),
            atomNumbers(atoms),
            atomAtomicNumbers(atoms),
            atomElements(atoms),
            atomCoordinates(atoms),
            atomCharges(atoms),
            atomVisibility(atoms)};

        carbohydrate::ResidueData residueData {
            residueNames(residues), residueTypes(residues), residueStringIds(residues), residueNumbers(residues)};

        carbohydrate::EdgeData edgeData {std::vector<std::string>(edgeCount(atomGraph), "")};

        return {atomData, residueData, edgeData, indices.indices, graph.atoms.source};
    }

    assembly::Graph createVisibleAssemblyGraph(const carbohydrate::CarbohydrateData& data)
    {
        graph::Database atomGraph = data.atomGraph;
        atomGraph.nodes.alive = data.atoms.visible;
        return createAssemblyGraph(data.indices, atomGraph);
    }
} // namespace gmml
