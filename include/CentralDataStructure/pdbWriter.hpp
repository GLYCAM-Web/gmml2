#ifndef INCLUDE_CENTRALDATASTRUCTURE_PDBWRITER_HPP
#define INCLUDE_CENTRALDATASTRUCTURE_PDBWRITER_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/fileType/pdb/pdbFileData.hpp"

#include <ostream>
#include <string>
#include <vector>

// Not all the options are available now, like e.g. writing molecule chain numbers, but I will add them as I need them
namespace gmml
{
    pdb::PdbFileAtomData toPdbFileAtomData(const std::vector<Atom*>& atoms, std::vector<std::string> recordNames);
    pdb::PdbFileData toPdbFileData(const GraphObjects& objects);

    void writeTrajectoryToPdb(
        std::ostream& stream, const pdb::PdbData& data, const std::vector<size_t>& selectedResidues);
    void WritePdb(std::ostream& stream, const GraphIndexData& graphData, const std::vector<std::string>& headerLines);
} // namespace gmml
#endif
