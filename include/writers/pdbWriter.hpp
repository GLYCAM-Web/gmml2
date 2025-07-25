#ifndef INCLUDE_WRITERS_PDBWRITER_HPP
#define INCLUDE_WRITERS_PDBWRITER_HPP

#include "include/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/readers/Pdb/pdbData.hpp"

#include <ostream>
#include <string>
#include <vector>

// Not all the options are available now, like e.g. writing molecule chain numbers, but I will add them as I need them
namespace gmml
{
    std::vector<bool> residueTER(const std::vector<ResidueType>& types);
    PdbFileAtomData toPdbFileAtomData(const std::vector<Atom*>& atoms, std::vector<std::string> recordNames);
    PdbFileData toPdbFileData(const GraphObjects& objects);

    void writeTrajectoryToPdb(
        std::ostream& stream, const pdb::PdbData& data, const std::vector<size_t>& selectedResidues);
    void WritePdb(std::ostream& stream, const GraphIndexData& graphData, const std::vector<std::string>& headerLines);
} // namespace gmml
#endif
