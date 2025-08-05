#ifndef INCLUDE_FILETYPE_PREP_PREPFILE_HPP
#define INCLUDE_FILETYPE_PREP_PREPFILE_HPP

#include "include/fileType/prep/prepAtom.hpp"
#include "include/fileType/prep/prepDataTypes.hpp"
#include "include/fileType/prep/prepResidue.hpp"

#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    namespace prep
    {
        PrepData readPrepFile(const std::string& prep_file);
        PrepData readPrepFile(const std::string& prep_file, const std::vector<std::string> queryNames);
        void setAtomConnectivities(PrepData& data);
        void generate3dStructures(PrepData& data);
        void write(const PrepData& data, const std::string& prep_file);
        void write(const PrepData& data, std::ostream& out_stream);
        std::string print(const PrepData& data);

        void readAllResidues(PrepData& data, std::istream& in_file);
        void readQueryResidues(PrepData& data, std::istream& in_file, const std::vector<std::string>& queryNames);
    } // namespace prep
} // namespace gmml

#endif
