#ifndef INCLUDE_READERS_PREP_PREPFILE_HPP
#define INCLUDE_READERS_PREP_PREPFILE_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/readers/Prep/prepAtom.hpp"
#include "include/readers/Prep/prepResidue.hpp"

#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    namespace prep
    {
        void readPrepFile(
            Molecule* molecule, std::vector<PrepResidueProperties>& properties, const std::string& prep_file);
        void readPrepFile(
            Molecule* molecule,
            std::vector<PrepResidueProperties>& properties,
            const std::string& prep_file,
            const std::vector<std::string> queryNames);
        void setAtomConnectivities(Molecule* molecule, const std::vector<PrepResidueProperties>& properties);
        void generate3dStructures(Molecule* molecule, std::vector<PrepResidueProperties>& properties);
        void write(
            Molecule* molecule, const std::vector<PrepResidueProperties>& properties, const std::string& prep_file);
        void write(Molecule* molecule, const std::vector<PrepResidueProperties>& properties, std::ostream& out_stream);
        std::string print(Molecule* molecule, const std::vector<PrepResidueProperties>& properties);

        void readAllResidues(Molecule* molecule, std::vector<PrepResidueProperties>& properties, std::istream& in_file);
        void readQueryResidues(
            Molecule* molecule,
            std::vector<PrepResidueProperties>& properties,
            std::istream& in_file,
            const std::vector<std::string>& queryNames);
    } // namespace prep
} // namespace gmml

#endif
