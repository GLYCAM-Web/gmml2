#ifndef INCLUDES_READERS_PREP_PREPFILE_HPP
#define INCLUDES_READERS_PREP_PREPFILE_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/readers/Prep/prepAtom.hpp"
#include "include/readers/Prep/prepResidue.hpp"

#include <istream>
#include <map>
#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    namespace prep
    {
        void readPrepFile(Molecule* molecule, const std::string& prep_file);
        void readPrepFile(Molecule* molecule, const std::string& prep_file, const std::vector<std::string> queryNames);
        void SetAtomConnectivities(Molecule* molecule);
        void Generate3dStructures(Molecule* molecule);
        void Write(Molecule* molecule, const std::string& prep_file);
        void Write(Molecule* molecule, std::ostream& out_stream);
        std::string Print(Molecule* molecule);

        void ReadAllResidues(Molecule* molecule, std::istream& in_file);
        void ReadQueryResidues(Molecule* molecule, std::istream& in_file, const std::vector<std::string>& queryNames);
    } // namespace prep
} // namespace gmml

#endif
