#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06Functions.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06ResidueNameGenerator.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Shapers/atomToCoordinateInterface.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageCreation.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp" // serializeAtomNumbers
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include <sstream>
#include <fstream>
#include <cctype>    // isDigit
#include <algorithm> //  std::erase, std::remove

using cdsCondensedSequence::Carbohydrate;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Carbohydrate::Carbohydrate(std::string inputSequence) : SequenceManipulator {inputSequence}
{
    this->setName("CONDENSEDSEQUENCE");
    this->ReorderSequence(); // So output is consistent regardless of input order e.g. Fuca1-2[Gala1-3]Glca vs
                             // Gala1-3[Fuca1-2]Glca. Same 3D structure.
    cdsParameters::ParameterManager parameterManager(this->GetGlycamNamesOfResidues());
    // prep::PrepFile glycamPrepFileSelect(prepFilePath, this->GetGlycamNamesOfResidues());
    for (auto& cdsResidue : this->getResidues())
    { // Move atoms from prep file into parsedResidues.
        if (cdsResidue->GetType() != cds::ResidueType::Deoxy)
        {
            ParsedResidue* parsedResidue = static_cast<ParsedResidue*>(cdsResidue);
            parameterManager.createAtomsForResidue(cdsResidue, this->GetGlycamResidueName(parsedResidue));
            if (parsedResidue->GetType() == cds::ResidueType::Derivative)
            { // Deal with adjusting charges for derivatives
                this->DerivativeChargeAdjustment(parsedResidue);
            }
        }
    }
    for (auto& cdsResidue : this->getResidues())
    { // Apply any deoxy
        if (cdsResidue->GetType() == cds::ResidueType::Deoxy)
        {
            this->ApplyDeoxy(static_cast<ParsedResidue*>(cdsResidue));
        }
    }
    // Have atom numbers go from 1 to number of atoms. Note this should be after deleting atoms due to deoxy
    this->SetIndexByConnectivity(); // For reporting residue index numbers to the user
    cds::serializeNumbers(this->getAtoms());
    // Set 3D structure
    this->DepthFirstSetConnectivityAndGeometry(this->GetTerminal()); // recurve start with terminal
    // Re-numbering is a hack as indices have global scope and two instances give too high numbers.
    unsigned int linkageIndex = 0;
    // Linkages should be Edges to avoid this as they already get renumbered above.
    for (auto& linkage : glycosidicLinkages_) // These will exist on the vector in order of edge connectivity set above.
    { // Greedy first means the atoms-to-move needs to be updated for every linkage:
        linkage.index = linkageIndex++;
        cds::determineAtomsThatMove(linkage.rotatableDihedrals);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Final carbohydrate overlap resolution starting.");
    this->ResolveOverlaps();
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

void Carbohydrate::replaceAglycone(cds::Residue* newAglycone)
{
    for (cds::ResidueLinkage& linkage : glycosidicLinkages_)
    {
        auto& linkageResidues = linkage.link.residues;
        if (linkageResidues.first == this->GetAglycone() || linkageResidues.second == this->GetAglycone())
        { // Old aglycone atoms are connected to reducing residue, are found during creation of rotatable dihedrals.
            this->cds::Molecule::deleteResidue(this->GetAglycone());
            newAglycone->addNeighbor(newAglycone->getName() + "-" + this->GetReducingResidue()->getName(),
                                     this->GetReducingResidue());
            cds::ResidueLink link = cds::findResidueLink({this->GetReducingResidue(), newAglycone});
            linkage               = cds::createResidueLinkage(link);
            cds::setDefaultShapeUsingMetadata(linkage.rotatableDihedrals);
            return;
        }
    }
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
// std::string fileOutputDirectory = "unspecified", std::string fileType = "PDB", std::string outputFileNaming =
// "structure"
void Carbohydrate::Generate3DStructureFiles(std::string fileOutputDirectory, std::string outputFileNaming)
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
        // Pdb file
        std::string completeFileName = PathAndFileName + ".pdb";
        std::ofstream outFileStream;
        outFileStream.open(completeFileName.c_str());
        this->WritePdb(outFileStream);
        outFileStream.close();
        // Off file
        completeFileName = PathAndFileName + ".off";
        outFileStream.open(completeFileName.c_str());
        this->WriteOff(outFileStream);
        outFileStream.close();
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
    unsigned long long int numberOfShapes = 1;
    for (auto& linkage : glycosidicLinkages_)
    {
        auto& dihedrals = linkage.rotatableDihedrals;
        numberOfShapes  *= (likelyShapesOnly ? cds::numberOfLikelyShapes(dihedrals) : cds::numberOfShapes(dihedrals));
    }
    return numberOfShapes;
}

cds::Residue* Carbohydrate::GetReducingResidue()
{ // Kindly dumb, but centralized stupidity. Tagging during construction would be better. Ano-ano linkages won't work,
  // but this is used by gp builder so ok.
    for (auto residue : this->getResidues())
    { // Return the first sugar residue that isn't an Aglycone
        if (residue->GetType() != ResidueType::Aglycone && residue->GetType() == ResidueType::Sugar)
        {
            return residue;
        }
    }
    if (this->GetResidueCount() > 1)
    {
        return this->getResidues().at(1);
    }
    std::string message = "Reducing residue requested for Carbohydrate with name " + this->getName() +
                          ", but it doesn't have more than 1 residue";
    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
    throw std::runtime_error(message);
}

cds::Residue* Carbohydrate::GetAglycone()
{ // Kindly dumb, but centralized stupidity. Tagging during construction would be better. Ano-ano linkages won't work,
  // but this is used by gp builder so ok.
    for (auto residue : this->getResidues())
    {
        if (residue->GetType() == ResidueType::Aglycone)
        {
            return residue;
        }
    }
    if (this->GetResidueCount() > 0)
    {
        return this->getResidues().front();
    }
    std::string message = "Aglycone residue requested for Carbohydrate with name " + this->getName() +
                          ", but it doesn't have even 1 residue";
    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
    throw std::runtime_error(message);
}

cds::Atom* Carbohydrate::GetAnomericAtom()
{
    return cdsSelections::guessAnomericAtomByForeignNeighbor(this->GetReducingResidue());
}

//////////////////////////////////////////////////////////
//                  PRIVATE FUNCTIONS                   //
//////////////////////////////////////////////////////////
void Carbohydrate::ApplyDeoxy(ParsedResidue* deoxyResidue)
{
    ParsedResidue* residueToBeDeoxified = deoxyResidue->GetParent();
    residueToBeDeoxified->MakeDeoxy(deoxyResidue->GetLink());
    this->deleteResidue(deoxyResidue); // Remove the deoxy derivative now.
    return;
}

void Carbohydrate::DerivativeChargeAdjustment(ParsedResidue* parsedResidue)
{
    std::string adjustAtomName = GlycamMetadata::GetAdjustmentAtom(parsedResidue->getName());
    adjustAtomName             += parsedResidue->GetLinkageName().substr(0, 1);

    cds::Atom* atomToAdjust = parsedResidue->GetParent()->FindAtom(adjustAtomName);
    atomToAdjust->setCharge(atomToAdjust->getCharge() + GlycamMetadata::GetAdjustmentCharge(parsedResidue->getName()));
    return;
}

void Carbohydrate::ConnectAndSetGeometry(cds::Residue* childResidue, cds::Residue* parentResidue)
{
    using cds::Atom;
    using cds::ResidueType;
    std::string linkageLabel = static_cast<ParsedResidue*>(childResidue)->GetLinkageName();
    // This is using the new Node<Residue> functionality and the old AtomNode
    // Now go figure out how which Atoms to bond to each other in the residues.
    // Rule: Can't ever have a child aglycone or a parent derivative.
    Atom* parentAtom         = nullptr;
    std::string childAtomName, parentAtomName;
    if (parentResidue->GetType() == ResidueType::Aglycone)
    {
        parentAtomName = GlycamMetadata::GetConnectionAtomForResidue(parentResidue->getName());
        parentAtom     = parentResidue->FindAtom(parentAtomName);
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
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        if (!isdigit(linkageLabel.substr(linkPosition).at(0)))
        {
            std::string message = "Could not convert the last linkage number to an integer: " + linkageLabel;
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        parentAtom =
            cdsSelections::getNonCarbonHeavyAtomNumbered(parentResidue->getAtoms(), linkageLabel.substr(linkPosition));
    }
    else
    {
        std::string message = "Error: parent residue: " + parentResidue->getName() +
                              " isn't either Aglycone or Sugar, and derivatives cannot be parents.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    if (parentAtom == nullptr)
    {
        std::string message = "Did not find connection atom in residue: " + parentResidue->getStringId();
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    // Now get child atom
    if (childResidue->GetType() == ResidueType::Derivative)
    {
        std::string glycamNameForResidue = this->GetGlycamResidueName(static_cast<ParsedResidue*>(childResidue));
        childAtomName                    = GlycamMetadata::GetConnectionAtomForResidue(glycamNameForResidue);
    }
    else if (childResidue->GetType() == ResidueType::Sugar)
    {
        std::string childLinkageNumber = linkageLabel.substr(1, 1);
        if (!isdigit(childLinkageNumber.at(0)))
        {
            std::string message = "Could not convert the first linkage number to an integer: " + childLinkageNumber;
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        childAtomName = "C" + childLinkageNumber;
    }
    else
    {
        std::string message = "Error: child residue: " + childResidue->getName() +
                              " is neither derivative or Sugar (aglycones cannot be children)";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    Atom* childAtom = childResidue->FindAtom(childAtomName);
    if (childAtom == nullptr)
    {
        std::string message =
            "Did not find child atom named " + childAtomName + " in child residue: " + childResidue->getStringId();
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    // Geometry
    cds::moveConnectedAtomsAccordingToBondLength(parentAtom, childAtom);
    //   Now bond the atoms. This could also set distance?, and angle? if passed to function?
    childAtom->addBond(parentAtom); // parentAtom also connected to childAtom. Fancy.
    for (auto& parentAtomNeighbor : parentAtom->getNeighbors())
    {
        if ((parentAtomNeighbor->getName().at(0) != 'H') && (parentAtomNeighbor != childAtom))
        {
            auto matrix =
                cds::rotationTo(std::array<Coordinate, 3> {*parentAtomNeighbor->getCoordinate(),
                                                           *parentAtom->getCoordinate(), *childAtom->getCoordinate()},
                                constants::degree2Radian(constants::DEFAULT_ANGLE));
            matrix.rotateCoordinates(childResidue->getCoordinates());
            break;
        }
    }
    // GREEDY: taken care of, but note that the atoms that move in RotatableDihedral class need to be updated after more
    // residues are added.
    cds::ResidueLink link        = cds::findResidueLink({childResidue, parentResidue});
    cds::ResidueLinkage& linkage = glycosidicLinkages_.emplace_back(cds::createResidueLinkage(link));
    cds::setDefaultShapeUsingMetadata(linkage.rotatableDihedrals);
    std::vector<cds::Atom*> childAtoms  = childResidue->getAtoms();  // keeps them alive in memory
    std::vector<cds::Atom*> parentAtoms = parentResidue->getAtoms(); // keeps them alive in memory
    cds::simpleWiggleCurrentRotamers(linkage.rotatableDihedrals, childAtoms, parentAtoms, 5);
    return;
}

// Gonna choke on cyclic glycans. Add a check for IsVisited when that is required.
void Carbohydrate::DepthFirstSetConnectivityAndGeometry(cds::Residue* currentParent)
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
        this->ConnectAndSetGeometry(child, currentParent);
        this->DepthFirstSetConnectivityAndGeometry(child);
    }
    return;
}

std::vector<std::string> Carbohydrate::GetGlycamNamesOfResidues() const
{
    std::vector<std::string> names(this->getResidues().size()); // set size of vec for speed.
    for (auto& residue : this->getResidues())
    {
        if (residue->GetType() != cds::ResidueType::Deoxy)
        {
            names.push_back(this->GetGlycamResidueName(static_cast<ParsedResidue*>(residue)));
        }
    }
    return names;
}

std::string Carbohydrate::GetGlycamResidueName(ParsedResidue* residue) const
{
    if (residue->GetType() == cds::ResidueType::Deoxy)
    {
        gmml::log(
            __LINE__, __FILE__, gmml::WAR,
            "Bad idea: We asked for Glycam Residue Name of a deoxy type residue (e.g. the 6D of Glc[6D]) with name: " +
                residue->GetResidueName());
        return "";
    }
    std::string linkages = "";
    if (residue->GetType() == cds::ResidueType::Sugar)
    {
        linkages = residue->GetChildLinkagesForGlycamResidueNaming();
    }
    std::string code = GlycamMetadata::Glycam06ResidueNameGenerator(
        linkages, residue->GetIsomer(), residue->GetResidueName(), residue->GetRingType(),
        residue->GetResidueModifier() + residue->GetRingShape(), residue->GetConfiguration());
    return code;
}

void Carbohydrate::SetDefaultShapeUsingMetadata()
{
    for (auto& linkage : glycosidicLinkages_)
    {
        cds::setDefaultShapeUsingMetadata(linkage.rotatableDihedrals);
    }
    return;
}

void Carbohydrate::ResolveOverlaps()
{
    for (auto& linkage : glycosidicLinkages_)
    {
        cds::simpleWiggleCurrentRotamers(linkage.rotatableDihedrals,
                                         {linkage.nonReducingOverlapResidues, linkage.reducingOverlapResidues}, 5);
    }
    return;
}
