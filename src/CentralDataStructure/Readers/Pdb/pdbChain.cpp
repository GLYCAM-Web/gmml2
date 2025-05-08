#include "includes/CentralDataStructure/Readers/Pdb/pdbChain.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/biology.hpp"

#include <string>
#include <vector>
#include <variant>

namespace
{
    void addResidueAtoms(pdb::PdbData& data, size_t residueId,
                         const std::vector<std::pair<std::string, cds::Coordinate>>& atoms)
    {
        for (auto& atom : atoms)
        {
            const std::string& name      = atom.first;
            const cds::Coordinate& coord = atom.second;
            pdb::addPdbAtom(data, residueId, name, coord);
        }
    }

    void addResidueBonds(pdb::PdbData& data, size_t residueId,
                         const std::vector<std::array<std::variant<size_t, std::string>, 2>>& bonds)
    {
        auto atomIndex = [&](const std::variant<size_t, std::string>& key)
        {
            if (std::holds_alternative<size_t>(key))
            {
                return std::get<size_t>(key);
            }
            else
            {
                return findResidueAtom(data, residueId, std::get<std::string>(key));
            }
        };
        for (auto& bond : bonds)
        {
            pdb::addBond(data, atomIndex(bond[0]), atomIndex(bond[1]));
        }
    }
} // namespace

void pdb::readChain(PdbData& data, size_t moleculeId, std::stringstream& stream_block)
{
    std::string line;
    while (getline(stream_block, line))
    {
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
        if ((recordName == "ATOM") || (recordName == "HETATM"))
        {
            std::stringstream singleResidueSection = extractSingleResidueFromRecordSection(stream_block, line);
            size_t residueId                       = data.indices.residueMolecule.size();
            data.indices.residueMolecule.push_back(moleculeId);
            data.moleculeResidueOrder[moleculeId].push_back(residueId);
            cds::Residue* residue = data.indices.molecules[moleculeId]->addResidue(std::make_unique<cds::Residue>());
            data.indices.residues.push_back(residue);
            readResidue(data, residueId, singleResidueSection, line);
        }
        else
        {
            gmml::log(__LINE__, __FILE__, gmml::WAR,
                      "In PdbChain Constructor with record that isn't cool: " + recordName);
            break;
        }
    }
    tagTerminalResidues(data, moleculeId);
}

void pdb::tagTerminalResidues(PdbData& data, size_t moleculeId)
{
    size_t nTer = getNTerminal(data, moleculeId);
    if (nTer < data.indices.residues.size())
    {
        data.indices.residues[nTer]->setNTerminal();
    }
    size_t cTer = getCTerminal(data, moleculeId);
    if (cTer < data.indices.residues.size())
    {
        data.indices.residues[cTer]->setCTerminal();
    }
}

std::stringstream pdb::extractSingleResidueFromRecordSection(std::stringstream& pdbFileStream, std::string line)
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
    return singleResidueSection;
}

void pdb::InsertCap(PdbData& data, size_t moleculeId, size_t refResidueId, const std::string& type)
{
    // This approach is bad. When parameter manager is good we can use that to remove the get_carestian stuff
    using cds::Coordinate;
    auto coord = [&](size_t n)
    {
        return data.indices.atoms[n]->coordinate();
    };
    cds::Molecule* molecule  = data.indices.molecules[moleculeId];
    cds::Residue* refResidue = data.indices.residues[refResidueId];
    if (type == "NHCH3") // NME
    {
        //        NME resid numbers. Otherwise good.
        Coordinate cCoordProtein  = coord(findResidueAtom(data, refResidueId, "C"));
        Coordinate caCoordProtein = coord(findResidueAtom(data, refResidueId, "CA"));
        Coordinate oCoordProtein  = coord(findResidueAtom(data, refResidueId, "O"));
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
        size_t residueId = data.indices.residueMolecule.size();
        size_t refIndex  = codeUtils::indexOf(data.moleculeResidueOrder[moleculeId], refResidueId);
        data.indices.residueMolecule.push_back(moleculeId);
        std::vector<size_t>& residueOrder = data.moleculeResidueOrder[moleculeId];
        residueOrder.insert(residueOrder.begin() + refIndex + 1, residueId);
        cds::Residue* newNMEResidue =
            molecule->insertNewResidue(std::make_unique<cds::Residue>("NME", refResidue), *refResidue);
        data.indices.residues.push_back(newNMEResidue);
        readResidue(data, residueId, refResidueId);
        addResidueAtoms(data, residueId,
                        {
                            {   "N",    nCoordNME},
                            {   "H",    hCoordNME},
                            { "CH3",  ch3CoordNME},
                            {"HH31", hh31CoordNME},
                            {"HH32", hh32CoordNME},
                            {"HH33", hh33CoordNME}
        });
        addResidueBonds(data, residueId,
                        {
                            {"N", findResidueAtom(data, refResidueId, "C")},
                            {"N", "H"},
                            {"N", "CH3"},
                            {"CH3", "HH31"},
                            {"CH3", "HH32"},
                            {"CH3", "HH33"}
        });
        newNMEResidue->SetType(cds::ResidueType::ProteinCappingGroup);
        data.residues.hasTerCard[residueId] = true;
    }
    else if (type == "COCH3") // ACE
    {
        //        NME resid numbers. Otherwise good.
        // These are the atoms in residue that I use to build the ACE out from.
        Coordinate cCoordProtein  = coord(findResidueAtom(data, refResidueId, "C"));
        Coordinate caCoordProtein = coord(findResidueAtom(data, refResidueId, "CA"));
        Coordinate nCoordProtein  = coord(findResidueAtom(data, refResidueId, "N"));
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
        auto refPosition = molecule->findPositionOfResidue(refResidue);
        --refPosition;
        cds::Residue* previousResidue =
            (*refPosition).get(); // Its an iterator to a unique ptr, so deref and get the raw. Ugh.
        size_t residueId = data.indices.residueMolecule.size();
        size_t refId     = codeUtils::indexOf(data.indices.residues, previousResidue);
        size_t refIndex  = codeUtils::indexOf(data.moleculeResidueOrder[moleculeId], refId);
        data.indices.residueMolecule.push_back(moleculeId);
        std::vector<size_t>& residueOrder = data.moleculeResidueOrder[moleculeId];
        residueOrder.insert(residueOrder.begin() + refIndex + 1, residueId);
        cds::Residue* newACEResidue =
            molecule->insertNewResidue(std::make_unique<cds::Residue>("ACE", previousResidue), *previousResidue);
        data.indices.residues.push_back(newACEResidue);
        readResidue(data, residueId, refId);
        addResidueAtoms(data, residueId,
                        {
                            {   "C",    cCoordACE},
                            {   "O",    oCoordACE},
                            { "CH3",  ch3CoordACE},
                            {"HH31", hh31CoordACE},
                            {"HH32", hh32CoordACE},
                            {"HH33", hh33CoordACE}
        });
        addResidueBonds(data, residueId,
                        {
                            {"C", findResidueAtom(data, refResidueId, "N")},
                            {"C", "O"},
                            {"C", "CH3"},
                            {"CH3", "HH31"},
                            {"CH3", "HH32"},
                            {"CH3", "HH33"}
        });
        newACEResidue->SetType(cds::ResidueType::ProteinCappingGroup);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Created ACE residue: " + pdbResidueId(data, residueId).print());
    }
}

void pdb::ModifyTerminal(PdbData& data, size_t residueId, const std::string& type)
{
    if (type == "NH3+") // For now, leaving it to tleap to add the correct H's
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying N Terminal of : " + pdbResidueId(data, residueId).print());
        size_t atomId = findResidueAtom(data, residueId, "H");
        deletePdbAtom(data, residueId, atomId);
    }
    else if (type == "CO2-")
    {
        size_t atomCount = data.indices.atoms.size();
        gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying C Terminal of : " + pdbResidueId(data, residueId).print());
        size_t atomOXT = findResidueAtom(data, residueId, "OXT");
        if (atomOXT < atomCount)
        {
            cds::Atom* atom = data.indices.atoms[atomOXT];
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "OXT atom already exists: " + atom->getName() + "_" + std::to_string(atom->getNumber()));
            return;
        }
        // I don't like this, but at least it's somewhat contained:
        size_t atomCA = findResidueAtom(data, residueId, "CA");
        size_t atomC  = findResidueAtom(data, residueId, "C");
        size_t atomO  = findResidueAtom(data, residueId, "O");
        if (atomCA == atomCount || atomC == atomCount || atomO == atomCount)
        {
            gmml::log(
                __LINE__, __FILE__, gmml::WAR,
                "Cterminal residue missing an atoms named CA, C or O, cannot create an OXT atom for this residue.");
            return;
        }
        auto coord = [&](size_t n)
        {
            return data.indices.atoms[n]->coordinate();
        };
        cds::Coordinate oxtCoord =
            cds::calculateCoordinateFromInternalCoords(coord(atomCA), coord(atomC), coord(atomO), 120.0, 180.0, 1.25);
        size_t oxtAtom = addPdbAtom(data, residueId, "OXT", oxtCoord);
        addBond(data, oxtAtom, atomC);
        cds::Atom* atomOptr = data.indices.atoms[atomO];
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Created new atom named OXT after " + atomOptr->getName() + "_" +
                      std::to_string(atomOptr->getNumber()));
        return;
    }
    gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    return;
}

// Only makes sense for proteins.
// Assumes vector is populated from N-terminal to C-terminal.
size_t pdb::getNTerminal(const PdbData& data, size_t moleculeId)
{
    std::function<bool(const size_t&)> isProtein = [&](size_t id)
    {
        return codeUtils::contains(biology::proteinResidueNames, data.indices.residues[id]->getName());
    };
    std::vector<size_t> proteinResidues = codeUtils::vectorFilter(isProtein, data.moleculeResidueOrder[moleculeId]);
    if (proteinResidues.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Looked for terminal residue of chain with protein residues.");
        return data.indices.residues.size();
    }
    return proteinResidues.front();
}

size_t pdb::getCTerminal(const PdbData& data, size_t moleculeId)
{
    std::function<bool(const size_t&)> isProtein = [&](size_t id)
    {
        return codeUtils::contains(biology::proteinResidueNames, data.indices.residues[id]->getName());
    };
    std::vector<size_t> proteinResidues = codeUtils::vectorFilter(isProtein, data.moleculeResidueOrder[moleculeId]);
    if (proteinResidues.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Looked for terminal residue of chain with protein residues.");
        return data.indices.residues.size();
    }
    return proteinResidues.back();
}
