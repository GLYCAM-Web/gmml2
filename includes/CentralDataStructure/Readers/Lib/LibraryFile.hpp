#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_LIB_LIBRARYFILE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_LIB_LIBRARYFILE_HPP

#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/Readers/Lib/LibraryResidue.hpp"
#include <istream>

namespace lib
{

    class LibraryFile : public cds::Molecule
    {
      public:
        LibraryFile(const std::string& filePath);
        LibraryFile(const std::string& filePath, const std::vector<std::string> queryNames);

      private:
        std::stringstream ExtractUnitSection(std::istream& inputFileStream, const std::string unitName);
        void ParseInFileStream(std::istream& inputFileStream);
        void ParseQueryResidues(std::istream& inputFileStream, const std::vector<std::string>& queryNames);
    };

} // namespace lib
#endif
