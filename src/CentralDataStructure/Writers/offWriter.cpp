#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>
#include <string>
#include <iomanip>

namespace
{
    std::vector<size_t> atomIndices(const std::vector<cds::Atom*>& atoms, const std::vector<cds::Atom*>& reindex)
    {
        std::vector<size_t> result;
        result.reserve(reindex.size());
        for (auto& find : reindex)
        {
            result.push_back(cds::atomVectorIndex(atoms, find));
        }
        return result;
    }

    std::vector<std::pair<size_t, size_t>> uniqueAtomBonds(const std::vector<cds::Atom*>& atoms)
    {
        std::vector<std::pair<size_t, size_t>> result;
        result.reserve(2 * atoms.size());
        for (auto& atom : atoms)
        {
            size_t index = cds::atomVectorIndex(atoms, atom);
            for (auto& neighbor : atom->getChildren())
            {
                size_t neighborIndex = cds::atomVectorIndex(atoms, neighbor);
                result.push_back({index, neighborIndex});
            }
        }
        return result;
    }
} // namespace

cds::OffFileData cds::toOffFileData(const std::vector<Residue*>& residues)
{
    std::vector<Atom*> atoms;
    std::vector<std::vector<size_t>> indices;
    std::vector<std::vector<size_t>> connections;
    std::vector<size_t> atomResidues;
    size_t residueIndex = 0;
    for (auto& residue : residues)
    {
        std::vector<Atom*> residueAtoms = residue->getAtoms();
        indices.push_back(codeUtils::indexVectorWithOffset(atoms.size(), residueAtoms));
        std::vector<Atom*> atomsConnected = atomsConnectedToOtherResidues(residueAtoms);
        codeUtils::insertInto(atoms, residueAtoms);
        codeUtils::insertInto(atomResidues, std::vector<size_t>(residueAtoms.size(), residueIndex));
        connections.push_back(atomIndices(atoms, atomsConnected));
        residueIndex++;
    }
    OffFileAtomData atomData {atomNumbers(atoms), atomNames(atoms),       atomTypes(atoms), atomAtomicNumbers(atoms),
                              atomCharges(atoms), atomCoordinates(atoms), atomResidues,     uniqueAtomBonds(atoms)};
    OffFileResidueData residueData {residueNumbers(residues), residueNames(residues), residueTypes(residues), indices,
                                    connections};
    return OffFileData {residueData, atomData};
}

void cds::serializeResiduesIndividually(std::vector<cds::Residue*>& residues)
{
    for (auto& residue : residues)
    {
        residue->setNumber(1);
        cds::serializeNumbers(residue->getAtoms());
    }
}

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

void cds::WriteOffFileUnit(const std::vector<size_t>& residueIndices, const OffFileResidueData& residues,
                           const OffFileAtomData& atoms, std::ostream& stream, const std::string& unitName)
{
    // WriteAtomSection
    const std::string FLAG = "131072";
    stream << "!entry." << unitName
           << ".unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg"
           << std::endl;
    for (size_t residueIndex : residueIndices)
    {
        unsigned int atomNumberInResidue = 1;
        for (size_t atomIndex : residues.atomIndices[residueIndex])
        {
            stream << " \"" << atoms.names[atomIndex] << "\" "
                   << "\"" << atoms.types[atomIndex] << "\" "
                   << "0"
                   << " " << residues.numbers[residueIndex] << " " << FLAG << " " << atomNumberInResidue << " "
                   << atoms.atomicNumbers[atomIndex] << " " << std::fixed << std::setprecision(6)
                   << atoms.charges[atomIndex] << std::endl;
            atomNumberInResidue++;
        }
    }
    // WriteAtomPertInfoSection
    stream << "!entry." << unitName
           << ".unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" << std::endl;
    for (size_t residueIndex : residueIndices)
    {
        for (size_t atomIndex : residues.atomIndices[residueIndex])
        {
            stream << " \"" << atoms.names[atomIndex] << "\" "
                   << "\"" << atoms.types[atomIndex] << "\" " << 0 << " " << -1 << " " << std::setprecision(1) << 0.0
                   << std::endl;
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
    stream << " " << residueIndices.size() + 1 << std::endl;
    // WriteConnectSection
    //  Note: this is silly but fine for most cases. If you're reading this it's because it mattered and you need to
    //  make it better.
    stream << "!entry." << unitName << ".unit.connect array int" << std::endl;
    stream << " " << 1 << std::endl;
    stream << " " << atoms.numbers[residues.atomIndices[residueIndices.back()].back()] << std::endl;
    // WriteConnectivitySection
    stream << "!entry." << unitName << ".unit.connectivity table  int atom1x  int atom2x  int flags" << std::endl;
    std::vector<bool> residueIncluded = codeUtils::indexMask(residues.names.size(), residueIndices);
    for (auto& bond : atoms.bonds)
    {
        if (residueIncluded[atoms.residues[bond.first]])
        {
            int number         = atoms.numbers[bond.first];
            int neighborNumber = atoms.numbers[bond.second];
            int min            = std::min(number, neighborNumber);
            int max            = std::max(number, neighborNumber);
            // According to docs: (the *second* atom is the one with the larger index). So ordering
            stream << " " << min << " " << max << " " << 1 << std::endl;
        }
    }
    // WriteHierarchySection
    stream << "!entry." << unitName << ".unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx"
           << std::endl;
    for (size_t residueIndex : residueIndices)
    {
        stream << " \""
               << "U"
               << "\""
               << " " << 0 << " "
               << "\""
               << "R"
               << "\""
               << " " << residues.numbers[residueIndex] << std::endl;
        for (size_t atomIndex : residues.atomIndices[residueIndex])
        {
            stream << " \""
                   << "R"
                   << "\""
                   << " " << residues.numbers[residueIndex] << " "
                   << "\""
                   << "A"
                   << "\""
                   << " " << atoms.numbers[atomIndex] << std::endl;
        }
    }
    // WriteNameSection
    stream << "!entry." << unitName << ".unit.name single str" << std::endl;
    stream << " \"" << unitName << "\"" << std::endl;
    // WritePositionSection
    stream << "!entry." << unitName << ".unit.positions table  dbl x  dbl y  dbl z" << std::endl;
    for (size_t residueIndex : residueIndices)
    {
        for (size_t atomIndex : residues.atomIndices[residueIndex])
        {
            Coordinate coord = atoms.coordinates[atomIndex];
            stream << std::setprecision(6) << std::fixed << " " << coord.GetX() << " " << coord.GetY() << " "
                   << coord.GetZ() << std::endl;
        }
    }
    // WriteResidueConnectSection // Every residue needs a head/tail regardless of reality. tleap uses this info.
    stream << "!entry." << unitName
           << ".unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x" << std::endl;
    for (size_t residueIndex : residueIndices)
    {
        std::vector<size_t> connectedAtoms = residues.atomsConnectedToOtherResidues[residueIndex];
        // Deal with residues that don't have a tail/head in reality:
        if (connectedAtoms.size() == 1)
        { // Repeating the same atom changes the tree structure in the parm7 file. Not sure anything uses that. Old gmml
          // code repeats so doing that.
            // For reducing terminal old code puts a 1 2 0 0 0 0. So not repeat. Changing to 2 2 0 0 0 0 causes first
            // atom not to be M (main). Might be an issue let's see.
            connectedAtoms.push_back(connectedAtoms.front());
        }
        for (auto& atomIndex : connectedAtoms)
        {
            stream << " " << atoms.numbers[atomIndex];
        }
        int columnsWithZero = 6 - connectedAtoms.size();
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
    for (size_t residueIndex : residueIndices)
    {
        const std::vector<size_t>& atomIndices = residues.atomIndices[residueIndex];
        unsigned int childseq                  = atomIndices.size() + 1;
        unsigned int startatomx                = atoms.numbers[atomIndices.front()];
        std::string restype                    = cds::getOffType(residues.types[residueIndex]);
        unsigned int imagingx                  = 0;
        stream << " \"" << residues.names[residueIndex] << "\""
               << " " << residues.numbers[residueIndex] << " " << childseq << " " << startatomx << " "
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
    for (size_t residueIndex : residueIndices)
    {
        for (size_t n = 0; n < residues.atomIndices[residueIndex].size(); n++)
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

void cds::WriteResiduesIndividuallyToOffFile(std::ostream& stream, const OffFileData& data)
{ // For writing each residue separately
    size_t residueCount = data.residues.names.size();
    stream << "!!index array str" << std::endl;
    for (size_t n = 0; n < residueCount; n++)
    {
        stream << " \"" << data.residues.names[n] << "\"" << std::endl;
    }
    for (size_t n = 0; n < residueCount; n++)
    {
        cds::WriteOffFileUnit({n}, data.residues, data.atoms, stream, data.residues.names[n]);
    }
}

void cds::WriteResiduesTogetherToOffFile(std::ostream& stream, const OffFileData& data, const std::string& unitName)
{ // For writing residues together as a molecule
    stream << "!!index array str" << std::endl;
    stream << " \"" << unitName << "\"" << std::endl;
    cds::WriteOffFileUnit(codeUtils::indexVector(data.residues.names), data.residues, data.atoms, stream, unitName);
}
