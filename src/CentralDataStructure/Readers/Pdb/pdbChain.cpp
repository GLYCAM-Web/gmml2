#include "includes/CentralDataStructure/Readers/Pdb/pdbChain.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp" // get_cartesian_point_from_internal_coords
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/biology.hpp" // proteinResidueNames

using pdb::PdbChain;

////////////////////////////////////////////////////////////
////                       CONSTRUCTOR                    //
////////////////////////////////////////////////////////////
PdbChain::PdbChain(std::stringstream& stream_block, const std::string& chainId) : cds::Molecule(chainId)
{
    // gmml::log(__LINE__, __FILE__, gmml::INF,
    //           "Constructing PdbChain from stream_block >>>>>>>>>:\n" + stream_block.str() +
    //               "\n<<<<<<<<<<<<<< end stream_block");
    std::string line;
    while (getline(stream_block, line))
    {
        // gmml::log(__LINE__,__FILE__,gmml::INF, "Constructing chain with line: " + line);
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
        if ((recordName == "ATOM") || (recordName == "HETATM"))
        {
            std::stringstream singleResidueSection = this->extractSingleResidueFromRecordSection(stream_block, line);
            this->addResidue(std::make_unique<PdbResidue>(singleResidueSection, line));
        }
        else
        {
            gmml::log(__LINE__, __FILE__, gmml::WAR,
                      "In PdbChain Constructor with record that isn't cool: " + recordName);
            break;
        }
    }
    // gmml::log(__LINE__, __FILE__, gmml::INF, "Adding NTerminal and CTerminal tags if protein present");
    this->tagTerminalResidues();
    // gmml::log(__LINE__, __FILE__, gmml::INF, "PdbChain Constructor Complete Captain");
    return;
}

////////////////////////////////////////////////////////////
////                         ACCESSOR                     //
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
////                    FUNCTIONS                         //
////////////////////////////////////////////////////////////

void PdbChain::tagTerminalResidues()
{
    PdbResidue* nTer = this->getNTerminal();
    if (nTer != nullptr)
    {
        nTer->addLabel("NTerminal");
    }
    PdbResidue* cTer = this->getCTerminal();
    if (cTer != nullptr)
    {
        cTer->addLabel("CTerminal");
    }
    return;
}

std::stringstream PdbChain::extractSingleResidueFromRecordSection(std::stringstream& pdbFileStream, std::string line)
{
    std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
    std::stringstream singleResidueSection;
    pdb::ResidueId residueId(line);
    pdb::ResidueId initialResidueId = residueId;
    while (residueId == initialResidueId)
    {
        singleResidueSection << line << std::endl;
        previousLinePosition = pdbFileStream.tellg(); // Save current line position.
        if (!std::getline(pdbFileStream, line))
        {
            break; // // If we hit the end, time to leave.
        }
        residueId = ResidueId(line);
    }
    // Go back to previous line position. E.g. was reading HEADER and found TITLE:
    pdbFileStream.seekg(previousLinePosition);
    //    gmml::log(__LINE__, __FILE__, gmml::INF,
    //              "Single residue section is:\n" + singleResidueSection.str() + "\nEnd of single residue section.");
    return singleResidueSection;
}

void PdbChain::InsertCap(const PdbResidue& refResidue, const std::string& type)
{
    // This approach is bad. When parameter manager is good we can use that to remove the get_carestian stuff
    using cds::Coordinate;
    if (type == "NHCH3") // NME
    {
        //        int sequenceNumber = refResidue.GetSequenceNumber() + 1; // Single gaps will end up with the same ACE
        //        NME resid numbers. Otherwise good.
        Coordinate cCoordProtein  = refResidue.FindAtom("C")->coordinate();
        Coordinate caCoordProtein = refResidue.FindAtom("CA")->coordinate();
        Coordinate oCoordProtein  = refResidue.FindAtom("O")->coordinate();
        Coordinate nCoordNME =
            cds::calculateCoordinateFromInternalCoords(oCoordProtein, caCoordProtein, cCoordProtein, 120.0, 180.0, 1.4);
        Coordinate hCoordNME =
            cds::calculateCoordinateFromInternalCoords(oCoordProtein, caCoordProtein, nCoordNME, 109.0, 180.0, 1.0);
        Coordinate ch3CoordNME =
            cds::calculateCoordinateFromInternalCoords(caCoordProtein, cCoordProtein, nCoordNME, 125.0, 180.0, 1.48);
        Coordinate hh31CoordNME =
            cds::calculateCoordinateFromInternalCoords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, 180.0, 1.09);
        Coordinate hh32CoordNME =
            cds::calculateCoordinateFromInternalCoords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, 60.0, 1.09);
        Coordinate hh33CoordNME =
            cds::calculateCoordinateFromInternalCoords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, -60.0, 1.09);
        cds::Residue* newNMEResidue =
            this->insertNewResidue(std::make_unique<PdbResidue>("NME", &refResidue), refResidue);
        cds::Atom* nAtom    = newNMEResidue->addAtom(std::make_unique<PdbAtom>("N", nCoordNME));
        cds::Atom* hAtom    = newNMEResidue->addAtom(std::make_unique<PdbAtom>("H", hCoordNME));
        cds::Atom* ch3Atom  = newNMEResidue->addAtom(std::make_unique<PdbAtom>("CH3", ch3CoordNME));
        cds::Atom* hh31Atom = newNMEResidue->addAtom(std::make_unique<PdbAtom>("HH31", hh31CoordNME));
        cds::Atom* hh32Atom = newNMEResidue->addAtom(std::make_unique<PdbAtom>("HH32", hh32CoordNME));
        cds::Atom* hh33Atom = newNMEResidue->addAtom(std::make_unique<PdbAtom>("HH33", hh33CoordNME));
        addBond(nAtom, refResidue.FindAtom("C"));
        addBond(nAtom, hAtom);
        addBond(nAtom, ch3Atom);
        addBond(ch3Atom, hh31Atom);
        addBond(ch3Atom, hh32Atom);
        addBond(ch3Atom, hh33Atom);
        newNMEResidue->SetType(cds::ResidueType::ProteinCappingGroup);
        codeUtils::erratic_cast<PdbResidue*>(newNMEResidue)->AddTerCard(); // No longer used?
    }
    else if (type == "COCH3") // ACE
    {
        //        int sequenceNumber = refResidue.GetSequenceNumber() - 1; // Single gaps will end up with the same ACE
        //        NME resid numbers. Otherwise good.
        // These are the atoms in residue that I use to build the ACE out from.
        Coordinate cCoordProtein  = refResidue.FindAtom("C")->coordinate();
        Coordinate caCoordProtein = refResidue.FindAtom("CA")->coordinate();
        Coordinate nCoordProtein  = refResidue.FindAtom("N")->coordinate();
        // This is bad, should use templates loaded from lib/prep file instead.
        Coordinate cCoordACE = cds::calculateCoordinateFromInternalCoords(cCoordProtein, caCoordProtein, nCoordProtein,
                                                                          120.0, -130.0, 1.4);
        Coordinate oCoordACE =
            cds::calculateCoordinateFromInternalCoords(caCoordProtein, nCoordProtein, cCoordACE, 120.0, 0.0, 1.23);
        Coordinate ch3CoordACE =
            cds::calculateCoordinateFromInternalCoords(caCoordProtein, nCoordProtein, cCoordACE, 125.0, 180.0, 1.48);
        Coordinate hh31CoordACE =
            cds::calculateCoordinateFromInternalCoords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 180.0, 1.09);
        Coordinate hh32CoordACE =
            cds::calculateCoordinateFromInternalCoords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 60.0, 1.09);
        Coordinate hh33CoordACE =
            cds::calculateCoordinateFromInternalCoords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, -60.0, 1.09);
        // Ok this next bit is convoluted, but I look up the position of the first atom in the protein residue and
        // insert the new Atom before it, and get passed back the position of the newly created atom, so I can use that
        // when creating the next one and so on. With ACE we want to insert before the residue, so I'm finding the
        // residue before here:
        auto refPosition = this->findPositionOfResidue(&refResidue);
        --refPosition;
        PdbResidue* previousResidue = codeUtils::erratic_cast<PdbResidue*>(
            (*refPosition).get()); // Its an iterator to a unique ptr, so deref and get the raw. Ugh.
        cds::Residue* newACEResidue =
            this->insertNewResidue(std::make_unique<PdbResidue>("ACE", previousResidue), *previousResidue);
        cds::Atom* cAtom    = newACEResidue->addAtom(std::make_unique<PdbAtom>("C", cCoordACE));
        cds::Atom* oAtom    = newACEResidue->addAtom(std::make_unique<PdbAtom>("O", oCoordACE));
        cds::Atom* ch3Atom  = newACEResidue->addAtom(std::make_unique<PdbAtom>("CH3", ch3CoordACE));
        cds::Atom* hh31Atom = newACEResidue->addAtom(std::make_unique<PdbAtom>("HH31", hh31CoordACE));
        cds::Atom* hh32Atom = newACEResidue->addAtom(std::make_unique<PdbAtom>("HH32", hh32CoordACE));
        cds::Atom* hh33Atom = newACEResidue->addAtom(std::make_unique<PdbAtom>("HH33", hh33CoordACE));
        addBond(cAtom, refResidue.FindAtom("N"));
        addBond(cAtom, oAtom);
        addBond(cAtom, ch3Atom);
        addBond(ch3Atom, hh31Atom);
        addBond(ch3Atom, hh32Atom);
        addBond(ch3Atom, hh33Atom);
        newACEResidue->SetType(cds::ResidueType::ProteinCappingGroup);
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Created ACE residue: " + codeUtils::erratic_cast<PdbResidue*>(newACEResidue)->printId());
    }
}

void PdbChain::ModifyTerminal(const std::string& type, PdbResidue* terminalResidue)
{
    if (type == "NH3+") // For now, leaving it to tleap to add the correct H's
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying N Terminal of : " + terminalResidue->printId());
        cds::Atom* found = terminalResidue->FindAtom("H");
        if (found != nullptr)
        {
            const PdbAtom* atom = codeUtils::erratic_cast<const PdbAtom*>(terminalResidue->FindAtom("H"));
            gmml::log(__LINE__, __FILE__, gmml::INF, "Deleting atom with id: " + atom->GetId());
            terminalResidue->deleteAtom(atom);
        }
        return;
    }
    else if (type == "CO2-")
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying C Terminal of : " + terminalResidue->printId());
        const cds::Atom* atom = terminalResidue->FindAtom("OXT");
        if (atom != nullptr)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "OXT atom already exists: " + codeUtils::erratic_cast<const PdbAtom*>(atom)->GetId());
            return;
        }
        // I don't like this, but at least it's somewhat contained:
        const cds::Atom* atomCA = terminalResidue->FindAtom("CA");
        cds::Atom* atomC        = terminalResidue->FindAtom("C");
        const cds::Atom* atomO  = terminalResidue->FindAtom("O");
        if (atomCA == nullptr || atomC == nullptr || atomO == nullptr)
        {
            gmml::log(
                __LINE__, __FILE__, gmml::WAR,
                "Cterminal residue missing an atoms named CA, C or O, cannot create an OXT atom for this residue.");
            return;
        }
        cds::Coordinate oxtCoord = cds::calculateCoordinateFromInternalCoords(atomCA->coordinate(), atomC->coordinate(),
                                                                              atomO->coordinate(), 120.0, 180.0, 1.25);
        cds::Atom* oxtAtom       = terminalResidue->addAtom(std::make_unique<PdbAtom>("OXT", oxtCoord));
        addBond(oxtAtom, atomC);
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Created new atom named OXT after " + codeUtils::erratic_cast<const PdbAtom*>(atomO)->GetId());
        return;
    }
    gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    return;
}

// Only makes sense for proteins.
// Assumes vector is populated from N-terminal to C-terminal.
pdb::PdbResidue* PdbChain::getNTerminal()
{
    std::vector<cds::Residue*> proteinResidues = this->getResidues(biology::proteinResidueNames);
    if (proteinResidues.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Looked for terminal residue of chain with protein residues.");
        return nullptr;
    }
    return codeUtils::erratic_cast<PdbResidue*>(proteinResidues.front());
}

pdb::PdbResidue* PdbChain::getCTerminal()
{
    std::vector<cds::Residue*> proteinResidues = this->getResidues(biology::proteinResidueNames);
    if (proteinResidues.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Looked for terminal residue of chain with protein residues.");
        return nullptr;
    }
    return codeUtils::erratic_cast<PdbResidue*>(proteinResidues.back());
}
