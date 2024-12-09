#include "includes/CodeUtils/arguments.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>

namespace codeUtils
{
    std::string programName(char* argv[])
    {
        std::vector<std::string> arguments(argv, argv + 1);
        return arguments[0];
    }

    Arguments readArguments(int argc, char* argv[], const std::vector<ArgDef>& defs)
    {
        std::vector<std::string> arguments(argv + 1, argv + argc);

        std::vector<std::string> unnamed;
        std::vector<std::string> optionNames;
        std::vector<int> optionIds;
        std::vector<std::string> optionValues;
        std::vector<bool> defFound(defs.size(), false);
        std::vector<std::string> defFoundStr(defs.size(), "");

        auto addDef = [&](size_t index, const std::string& value, const std::string name)
        {
            if (defFound[index])
            {
                throw std::runtime_error("duplicate arguments: " + defFoundStr[index] + ", " + name);
            }
            defFound[index]    = true;
            defFoundStr[index] = name;
            optionNames.push_back(name);
            optionIds.push_back(defs[index].id);
            optionValues.push_back(value);
        };

        for (size_t n = 0; n < arguments.size(); n++)
        {
            const std::string& arg = arguments[n];
            if (startsWith(arg, "--"))
            {
                std::string name = arg.substr(2, arg.size() - 2);
                auto hasName     = [&](const ArgDef& def)
                {
                    return def.longName == name;
                };
                auto it = std::find_if(defs.begin(), defs.end(), hasName);
                if (it == defs.end())
                {
                    throw std::runtime_error("unknown argument " + arg);
                }
                size_t index      = it - defs.begin();
                const ArgDef& def = defs[index];
                switch (def.type)
                {
                    case ArgType::flag:
                        {
                            addDef(index, "", arg);
                            break;
                        }
                    case ArgType::option:
                        {
                            n++;
                            if (n >= arguments.size())
                            {
                                throw std::runtime_error("no value provided for argument " + arg);
                            }
                            addDef(index, arguments[n], arg);
                            break;
                        }
                    default:
                        throw std::runtime_error("erroneous argument type in " + arg + ". This is not a user error");
                }
            }
            else if (startsWith(arg, "-"))
            {
                std::string chars = arg.substr(1, arg.size() - 1);
                for (size_t k = 0; k < chars.size(); k++)
                {
                    char c = chars[k];
                    std::string dashArg {'-', c};
                    auto hasName = [&](const ArgDef& def)
                    {
                        return def.shortName == c;
                    };
                    auto it = std::find_if(defs.begin(), defs.end(), hasName);
                    if (it == defs.end())
                    {
                        throw std::runtime_error("unknown argument " + dashArg);
                    }
                    size_t index      = it - defs.begin();
                    const ArgDef& def = defs[index];
                    switch (def.type)
                    {
                        case ArgType::flag:
                            {
                                addDef(index, "", dashArg);
                                break;
                            }
                        case ArgType::option:
                            {
                                if (chars.size() > 1)
                                {
                                    throw std::runtime_error("cannot group non-flag argument " + dashArg + " with -" +
                                                             chars + "\ntry " + dashArg + " <" + def.nameOfValue + ">");
                                }
                                n++;
                                if (n >= arguments.size())
                                {
                                    throw std::runtime_error("no value provided for argument " + dashArg);
                                }
                                addDef(index, arguments[n], dashArg);
                                break;
                            }
                        default:
                            throw std::runtime_error("erroneous argument type in " + dashArg +
                                                     ". This is not a user error");
                    }
                }
            }
            else
            {
                unnamed.push_back(arg);
            }
        }
        return {unnamed, optionNames, optionIds, optionValues};
    }

    void validateArgumentCount(const Arguments& arguments, const std::vector<ArgDef>& defs)
    {
        size_t argCount = arguments.unnamed.size();
        size_t argMin   = 0;
        size_t argMax   = 0;
        for (auto& def : defs)
        {
            if (def.type == ArgType::unnamed)
            {
                argMin += def.requirement == ArgReq::required;
                argMax++;
            }
        }
        if (argCount < argMin || argCount > argMax)
        {
            std::string numStr = (argMin == argMax)
                                     ? std::to_string(argMin)
                                     : ("between " + std::to_string(argMin) + " and " + std::to_string(argMax));
            throw std::runtime_error("expected " + numStr + " unnamed arguments, got " +
                                     std::to_string(arguments.unnamed.size()));
        }
    }

    std::string helpString(const std::string& programName, const std::vector<ArgDef>& defs)
    {
        std::pair<std::string, std::string> emptyBrace {"", ""};
        std::pair<std::string, std::string> brackets {"[", "]"};
        auto braceType = [&](ArgReq req)
        {
            return (req == ArgReq::optional) ? brackets : emptyBrace;
        };
        std::string initial = "usage: " + programName + " ";
        std::string indent(initial.size(), ' ');
        size_t widthLimit = 60;
        std::ostringstream ss;
        ss << initial;
        size_t accum = 0;
        for (auto& def : defs)
        {
            if (def.type != ArgType::unnamed)
            {
                if (accum >= widthLimit)
                {
                    ss << "\n" << indent;
                    accum = 0;
                }
                bool hasShortName     = def.shortName != ' ';
                bool hasLongName      = def.longName != "";
                std::string valueStr  = (def.type == ArgType::option ? (" <" + def.nameOfValue + ">") : "");
                std::string shortStr  = hasShortName ? "-" + std::string {def.shortName} + valueStr : "";
                std::string longStr   = hasLongName ? "--" + def.longName + valueStr : "";
                std::string separator = (hasShortName && hasLongName ? " | " : "");
                std::pair<std::string, std::string> brace = braceType(def.requirement);
                std::string str                           = brace.first + shortStr + separator + longStr + brace.second;
                ss << str << " ";
                accum += str.size() + 1;
            }
        }
        ss << "\n" << indent;
        for (auto& def : defs)
        {
            if (def.type == ArgType::unnamed)
            {
                std::pair<std::string, std::string> brace = braceType(def.requirement);
                std::string str                           = brace.first + def.nameOfValue + brace.second;
                ss << str << " ";
            }
        }
        ss << "\n";
        return ss.str();
    }
} // namespace codeUtils
