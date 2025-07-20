#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP

#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/cdsTypes.hpp"

#include <ostream>
#include <string>
#include <vector>

// Not all the options are available now, like e.g. writing molecule chain numbers, but I will add them as I need them
namespace cds
{
    std::vector<bool> residueTER(const std::vector<ResidueType>& types);
    PdbFileAtomData toPdbFileAtomData(const std::vector<cds::Atom*>& atoms, std::vector<std::string> recordNames);
    PdbFileData toPdbFileData(const cds::GraphObjects& objects);

    void writeTrajectoryToPdb(
        std::ostream& stream, const pdb::PdbData& data, const std::vector<size_t>& selectedResidues);
    void WritePdb(std::ostream& stream, const GraphIndexData& graphData, const std::vector<std::string>& headerLines);
} // namespace cds
#endif
