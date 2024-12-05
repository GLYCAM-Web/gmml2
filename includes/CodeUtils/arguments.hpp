#ifndef INCLUDES_CODEUTILS_ARGUMENTS_HPP
#define INCLUDES_CODEUTILS_ARGUMENTS_HPP

#include <string>
#include <vector>

namespace codeUtils
{
    struct Arguments
    {
        std::vector<std::string> unnamed;
        std::vector<std::string> flags;
        std::vector<std::pair<std::string, std::string>> options;
    };

    Arguments readArguments(int argc, char* argv[]);
} // namespace codeUtils
#endif
