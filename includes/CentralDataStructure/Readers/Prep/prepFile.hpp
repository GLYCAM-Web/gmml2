#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PREP_PREPFILE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PREP_PREPFILE_HPP

#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepAtom.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepResidue.hpp"
#include <map>
#include <string>
#include <istream>
#include <ostream>
#include <vector>

namespace prep
{
    void readPrepFile(cds::Molecule* molecule, const std::string& prep_file);
    void readPrepFile(cds::Molecule* molecule, const std::string& prep_file, const std::vector<std::string> queryNames);
    void SetAtomConnectivities(cds::Molecule* molecule);
    void Generate3dStructures(cds::Molecule* molecule);
    void Write(cds::Molecule* molecule, const std::string& prep_file);
    void Write(cds::Molecule* molecule, std::ostream& out_stream);
    std::string Print(cds::Molecule* molecule);

    void ReadAllResidues(cds::Molecule* molecule, std::istream& in_file);
    void ReadQueryResidues(cds::Molecule* molecule, std::istream& in_file, const std::vector<std::string>& queryNames);
} // namespace prep
#endif
