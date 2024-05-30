#ifndef GMML_INCLUDES_CODEUTILS_FIND_HPP
#define GMML_INCLUDES_CODEUTILS_FIND_HPP

#include <string>
#include <unordered_map>

namespace codeUtils
{
    std::string FindStringInStringMap(const std::string& s, const std::unordered_map<std::string, std::string>& sMap);

} // namespace codeUtils
#endif // GMML_INCLUDES_CODEUTILS_STRINGS_HPP
