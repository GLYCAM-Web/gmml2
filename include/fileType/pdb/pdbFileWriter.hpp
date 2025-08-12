#ifndef INCLUDE_FILETYPE_PDB_PDBFILEWRITER_HPP
#define INCLUDE_FILETYPE_PDB_PDBFILEWRITER_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/fileType/pdb/pdbFileData.hpp"
#include "include/metadata/residueTypes.hpp"

#include <array>
#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        std::vector<bool> residueTER(const std::vector<ResidueType>& types);

        void writeAssemblyToPdb(
            std::ostream& stream,
            const PdbFileData& data,
            const assembly::Graph& graph,
            const std::vector<std::vector<size_t>>& residueIndices,
            const std::vector<std::vector<bool>>& residueTER,
            const std::vector<std::array<size_t, 2>>& connectionIndices);

        void writeMoleculeToPdb(
            std::ostream& stream,
            const PdbFileData& data,
            const assembly::Graph& graph,
            const std::vector<size_t>& residueIndices,
            const std::vector<bool>& residueTER);

        void writeTrajectoryToPdb(
            std::ostream& stream,
            const PdbData& data,
            const assembly::Graph& graph,
            const std::vector<size_t>& selectedResidues);

        void writeAtomToPdb(
            std::ostream& stream,
            const PdbFileFormat& format,
            const PdbFileResidueData& residues,
            size_t residueIndex,
            const PdbFileAtomData& atoms,
            size_t atomIndex);

        void writeConectCards(
            std::ostream& stream,
            const std::vector<uint>& atomNumbers,
            const std::vector<std::array<size_t, 2>>& connectionIndices);

        void theEnd(std::ostream& stream);
    } // namespace pdb
} // namespace gmml

#endif
