#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_LIB_LIBRARYFILE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_LIB_LIBRARYFILE_HPP

#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/Readers/Lib/LibraryResidue.hpp"
#include <sstream>

namespace lib
{

    class LibraryFile : public cds::Molecule
    {
      public:
        LibraryFile();
        LibraryFile(const std::string& filePath);
        LibraryFile(const std::string& filePath, const std::vector<std::string> queryNames);

      private:
        std::stringstream ExtractUnitSection(std::ifstream& inputFileStream, const std::string unitName);
        void ParseInFileStream(std::ifstream& inputFileStream);
        void ParseQueryResidues(std::ifstream& inputFileStream, const std::vector<std::string>& queryNames);
    };

} // namespace lib
#endif
