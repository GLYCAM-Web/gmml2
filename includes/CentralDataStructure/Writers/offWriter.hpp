#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include <vector>
#include <string>
#include <ostream>

namespace cds
{
    std::string getOffType(const cds::ResidueType queryType);
    void WriteOffFileUnit(std::vector<cds::Residue*> residues, std::ostream& stream, const std::string unitName);
    void WriteResiduesToOffFile(std::vector<cds::Residue*> residues, std::ostream& stream);
    void WriteMoleculeToOffFile(cds::Molecule* molecule, std::ostream& stream, const std::string unitName);
    void WriteAssemblyToOffFile(cds::Assembly* assembly, std::ostream& stream, const std::string unitName);
} // namespace cds
#endif
