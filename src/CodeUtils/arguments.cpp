#include "includes/CodeUtils/arguments.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <string>
#include <vector>
#include <stdexcept>

namespace codeUtils
{
    Arguments readArguments(int argc, char* argv[])
    {
        std::vector<std::string> arguments(argv + 1, argv + argc);

        std::vector<std::string> unnamed;
        std::vector<std::string> flags;
        std::vector<std::pair<std::string, std::string>> options;
        std::vector<std::string> errors;

        for (auto& arg : arguments)
        {
            if (startsWith(arg, "--"))
            {
                std::vector<std::string> nameValue = split(arg, '=');
                if (nameValue.size() == 1)
                {
                    flags.push_back(nameValue[0]);
                }
                else if (nameValue.size() == 2)
                {
                    options.push_back({nameValue[0], nameValue[1]});
                }
                else
                {
                    throw std::runtime_error("could not parse argument " + arg + "\n" +
                                             "expected --name or --name=value");
                }
            }
            else if (startsWith(arg, "-"))
            {
                throw std::runtime_error("could not parse argument " + arg + "\n" + "did you mean -" + arg);
            }
            else
            {
                unnamed.push_back(arg);
            }
        }
        return {unnamed, flags, options};
    }
} // namespace codeUtils
