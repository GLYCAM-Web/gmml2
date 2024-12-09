#ifndef INCLUDES_CODEUTILS_ARGUMENTS_HPP
#define INCLUDES_CODEUTILS_ARGUMENTS_HPP

#include <string>
#include <vector>

namespace codeUtils
{
    enum class ArgType
    {
        unnamed,
        flag,
        option
    };
    enum class ArgReq
    {
        optional,
        required
    };

    struct ArgDef
    {
        ArgReq requirement;
        ArgType type;
        int id;
        std::string longName    = "";
        char shortName          = ' ';
        std::string nameOfValue = "";
        std::string help        = "";
    };

    struct Arguments
    {
        std::vector<std::string> unnamed;
        std::vector<std::string> names;
        std::vector<int> ids;
        std::vector<std::string> values;
    };

    std::string programName(char* argv[]);
    Arguments readArguments(int argc, char* argv[], const std::vector<ArgDef>& defs);
    void validateArgumentCount(const Arguments& arguments, const std::vector<ArgDef>& defs);
    std::string helpString(const std::string& programName, const std::vector<ArgDef>& defs);
} // namespace codeUtils
#endif
