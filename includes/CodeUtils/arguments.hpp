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
        std::string name         = "";
        std::string defaultValue = "";
        std::string help         = "";
    };

    struct Arguments
    {
        std::string programName;
        std::vector<std::string> unnamed;
        std::vector<std::string> names;
        std::vector<std::string> values;
        std::vector<bool> hasValue;
    };

    Arguments readArguments(int argc, char* argv[]);
    void validateArgumentCount(const Arguments& arguments, const std::vector<ArgDef>& defs);
    void validateFlagsAndOptions(const Arguments& arguments, const std::vector<ArgDef>& defs);
    std::string helpString(const std::string& programName, const std::vector<ArgDef>& defs);
} // namespace codeUtils
#endif
