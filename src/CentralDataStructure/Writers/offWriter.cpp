#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <vector>
#include <string>
#include <iomanip>

namespace
{
    cds::AtomOffData toAtomOffData(const cds::Atom& atom)
    {
        std::vector<int> children = cds::atomNumbers(atom.getChildren());
        return cds::AtomOffData(atom.getNumber(), atom.getName(), atom.getType(), atom.getAtomicNumber(),
                                atom.getCharge(), atom.coordinate(), children);
    }

    cds::ResidueOffData toResidueOffData(const cds::Residue& residue)
    {
        std::vector<cds::Atom*> atoms = residue.getAtoms();
        std::vector<cds::AtomOffData> offAtoms;
        for (auto& atom : atoms)
        {
            offAtoms.push_back(toAtomOffData(*atom));
        }
        std::vector<cds::Atom*> atomsConnectedToOtherResidues = cds::atomsConnectedToOtherResidues(atoms);
        std::vector<int> connectedAtomNumbers                 = cds::atomNumbers(atomsConnectedToOtherResidues);
        return cds::ResidueOffData(residue.getNumber(), residue.getName(), residue.GetType(), offAtoms,
                                   connectedAtomNumbers);
    }

    std::vector<cds::ResidueOffData> toResidueOffDataVector(const std::vector<cds::Residue*>& residues)
    {
        std::vector<cds::ResidueOffData> result;
        for (auto& residue : residues)
        {
            result.push_back(toResidueOffData(*residue));
        }
        return result;
    }
} // namespace

cds::AtomOffData::AtomOffData(int number_, std::string name_, std::string type_, int atomicNumber_, double charge_,
                              Coordinate coordinate_, std::vector<int> children_)
    : number(number_), name(name_), type(type_), atomicNumber(atomicNumber_), charge(charge_), coordinate(coordinate_),
      children(children_)
{}

cds::ResidueOffData::ResidueOffData(int number_, std::string name_, ResidueType type_, std::vector<AtomOffData> atoms_,
                                    std::vector<int> connections_)
    : number(number_), name(name_), type(type_), atoms(atoms_), atomsConnectedToOtherResidues(connections_)
{}

std::string cds::getOffType(const cds::ResidueType queryType)
{
    if (queryType == cds::ResidueType::Protein)
    {
        return "p";
    }
    if (queryType == cds::ResidueType::Solvent)
    {
        return "w";
    }
    return "?";
}

void cds::WriteOffFileUnit(const std::vector<ResidueOffData>& residues, std::ostream& stream,
                           const std::string unitName)
{
    // WriteAtomSection
    const std::string FLAG = "131072";
    stream << "!entry." << unitName
           << ".unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg"
           << std::endl;
    for (auto& residue : residues)
    {
        unsigned int atomNumberInResidue = 1;
        for (auto& atom : residue.atoms)
        {
            stream << " \"" << atom.name << "\" "
                   << "\"" << atom.type << "\" "
                   << "0"
                   << " " << residue.number << " " << FLAG << " " << atomNumberInResidue << " " << atom.atomicNumber
                   << " " << std::fixed << std::setprecision(6) << atom.charge << std::endl;
            atomNumberInResidue++;
        }
    }
    // WriteAtomPertInfoSection
    stream << "!entry." << unitName
           << ".unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" << std::endl;
    for (auto& residue : residues)
    {
        for (auto& atom : residue.atoms)
        {
            stream << " \"" << atom.name << "\" "
                   << "\"" << atom.type << "\" " << 0 << " " << -1 << " " << std::setprecision(1) << 0.0 << std::endl;
        }
    }
    // WriteBoundBoxSection
    stream << "!entry." << unitName << ".unit.boundbox array dbl" << std::endl;
    stream << " "
           << "-1.000000" << std::endl;
    stream << " " << 0.0 << std::endl;
    stream << " " << 0.0 << std::endl;
    stream << " " << 0.0 << std::endl;
    stream << " " << 0.0 << std::endl;
    // WriteChildSequenceSection
    stream << "!entry." << unitName << ".unit.childsequence single int" << std::endl;
    stream << " " << residues.size() + 1 << std::endl;
    // WriteConnectSection
    //  Note: this is silly but fine for most cases. If you're reading this it's because it mattered and you need to
    //  make it better.
    stream << "!entry." << unitName << ".unit.connect array int" << std::endl;
    stream << " " << 1 << std::endl;
    stream << " " << residues.back().atoms.back().number << std::endl;
    // WriteConnectivitySection
    stream << "!entry." << unitName << ".unit.connectivity table  int atom1x  int atom2x  int flags" << std::endl;
    for (auto& residue : residues)
    {
        for (auto& atom : residue.atoms)
        {
            for (auto& neighborNumber : atom.children)
            { // According to docs: (the *second* atom is the one with the larger index). So ordering
                std::pair<int, int> numbers = (atom.number < neighborNumber)
                                                  ? std::pair<int, int> {atom.number, neighborNumber}
                                                  : std::pair<int, int> {neighborNumber, atom.number};
                stream << " " << numbers.first << " " << numbers.second << " " << 1 << std::endl;
            }
        }
    }
    // WriteHierarchySection
    stream << "!entry." << unitName << ".unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx"
           << std::endl;
    for (auto& residue : residues)
    {
        stream << " \""
               << "U"
               << "\""
               << " " << 0 << " "
               << "\""
               << "R"
               << "\""
               << " " << residue.number << std::endl;
        for (auto& atom : residue.atoms)
        {
            stream << " \""
                   << "R"
                   << "\""
                   << " " << residue.number << " "
                   << "\""
                   << "A"
                   << "\""
                   << " " << atom.number << std::endl;
        }
    }
    // WriteNameSection
    stream << "!entry." << unitName << ".unit.name single str" << std::endl;
    stream << " \"" << unitName << "\"" << std::endl;
    // WritePositionSection
    stream << "!entry." << unitName << ".unit.positions table  dbl x  dbl y  dbl z" << std::endl;
    for (auto& residue : residues)
    {
        for (auto& atom : residue.atoms)
        {
            Coordinate coord = atom.coordinate;
            stream << std::setprecision(6) << std::fixed << " " << coord.GetX() << " " << coord.GetY() << " "
                   << coord.GetZ() << std::endl;
        }
    }
    // WriteResidueConnectSection // Every residue needs a head/tail regardless of reality. tleap uses this info.
    stream << "!entry." << unitName
           << ".unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x" << std::endl;
    for (auto& residue : residues)
    {
        std::vector<int> atomNumbers = residue.atomsConnectedToOtherResidues;
        // Deal with residues that don't have a tail/head in reality:
        if (atomNumbers.size() == 1)
        { // Repeating the same atom changes the tree structure in the parm7 file. Not sure anything uses that. Old gmml
          // code repeats so doing that.
            // For reducing terminal old code puts a 1 2 0 0 0 0. So not repeat. Changing to 2 2 0 0 0 0 causes first
            // atom not to be M (main). Might be an issue let's see.
            atomNumbers.push_back(atomNumbers.front());
        }
        for (auto& number : atomNumbers)
        {
            stream << " " << number;
        }
        int columnsWithZero = 6 - atomNumbers.size();
        for (int i = 0; i < columnsWithZero; ++i)
        {
            stream << " "
                   << "0";
        }
        stream << std::endl;
    }
    // WriteResiduesSection
    stream << "!entry." << unitName
           << ".unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx"
           << std::endl;
    for (auto& residue : residues)
    {
        unsigned int childseq   = residue.atoms.size() + 1;
        unsigned int startatomx = residue.atoms.front().number;
        std::string restype     = cds::getOffType(residue.type);
        unsigned int imagingx   = 0;
        stream << " \"" << residue.name << "\""
               << " " << residue.number << " " << childseq << " " << startatomx << " "
               << "\"" << restype << "\""
               << " " << imagingx << std::endl;
    }
    // WriteSolventCapSection
    stream << "!entry." << unitName << ".unit.solventcap array dbl" << std::endl;
    stream << " "
           << "-1.000000" << std::endl;
    stream << " "
           << "0.0" << std::endl;
    stream << " "
           << "0.0" << std::endl;
    stream << " "
           << "0.0" << std::endl;
    stream << " "
           << "0.0" << std::endl;
    // WriteVelocitiesSection
    stream << "!entry." << unitName << ".unit.velocities table  dbl x  dbl y  dbl z" << std::endl;
    for (auto& residue : residues)
    {
        for (size_t n = 0; n < residue.atoms.size(); n++)
        { // Maybe later we'll want to deal with atom velocities...
            stream << " "
                   << "0.0"
                   << " "
                   << "0.0"
                   << " "
                   << "0.0" << std::endl;
        }
    }
    return;
}

void cds::WriteResiduesToOffFile(std::vector<cds::Residue*> residues, std::ostream& stream)
{ // For writing each residue separately
    stream << "!!index array str" << std::endl;
    for (auto& residue : residues)
    {
        stream << " \"" << residue->getName() << "\"" << std::endl;
    }
    for (auto& residue : residues)
    {
        std::vector<Residue*> residues = {residue};
        std::vector<Atom*> atoms       = residue->getAtoms();
        cds::serializeNumbers(atoms);
        cds::serializeNumbers(residues);
        cds::WriteOffFileUnit(toResidueOffDataVector(residues), stream, residue->getName());
    }
}

void cds::WriteToOffFile(const std::vector<Residue*>& residues, std::ostream& stream, const std::string unitName)
{ // For writing residues together as a molecule
    stream << "!!index array str" << std::endl;
    stream << " \"" << unitName << "\"" << std::endl;
    cds::WriteOffFileUnit(toResidueOffDataVector(residues), stream, unitName);
}
