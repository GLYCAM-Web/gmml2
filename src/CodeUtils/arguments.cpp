#include "includes/CodeUtils/arguments.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>

namespace codeUtils
{
    Arguments readArguments(int argc, char* argv[])
    {
        std::vector<std::string> arguments(argv, argv + argc);
        std::string programName = arguments[0];

        std::vector<std::string> unnamed;
        std::vector<std::string> optionNames;
        std::vector<std::string> optionValues;
        std::vector<bool> optionHasValue;
        std::vector<std::string> errors;

        for (size_t n = 1; n < arguments.size(); n++)
        {
            const std::string& arg = arguments[n];
            if (startsWith(arg, "--"))
            {
                std::vector<std::string> nameValue = split(arg, '=');
                if (nameValue.size() == 1)
                {
                    optionNames.push_back(nameValue[0]);
                    optionValues.push_back("");
                    optionHasValue.push_back(false);
                }
                else if (nameValue.size() == 2)
                {
                    optionNames.push_back(nameValue[0]);
                    optionValues.push_back(nameValue[1]);
                    optionHasValue.push_back(true);
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
        return {programName, unnamed, optionNames, optionValues, optionHasValue};
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

    void validateFlagsAndOptions(const Arguments& arguments, const std::vector<ArgDef>& defs)
    {
        const std::vector<std::string>& names = arguments.names;
        std::vector<bool> nameFound(names.size(), false);
        for (auto& def : defs)
        {
            if (def.type != ArgType::unnamed)
            {
                bool required = def.requirement == ArgReq::required;
                auto iter     = std::find(names.begin(), names.end(), def.name);
                bool found    = iter != names.end();
                if (required && !found)
                {
                    throw std::runtime_error("missing required argument " + def.name);
                }
                else if (found)
                {
                    size_t index = iter - names.begin();
                    if (nameFound[index])
                    {
                        throw std::runtime_error("duplicate argument " + def.name);
                    }
                    bool hasValue = def.type == ArgType::option;
                    if (arguments.hasValue[index] != hasValue)
                    {
                        if (hasValue)
                        {
                            throw std::runtime_error("value required for argument " + def.name);
                        }
                        else
                        {
                            throw std::runtime_error("no value expected for argument " + def.name);
                        }
                    }
                    nameFound[index] = true;
                }
            }
        }
        for (size_t n = 0; n < nameFound.size(); n++)
        {
            if (!nameFound[n])
            {
                throw std::runtime_error("unknown argument " + names[n]);
            }
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
                const std::pair<std::string, std::string>& brace = braceType(def.requirement);
                std::string str                                  = brace.first + def.name +
                                  (def.type == ArgType::option ? ("=<" + def.nameOfValue + ">") : "") + brace.second;
                ss << str << " ";
                accum += str.size() + 1;
                if (accum >= widthLimit)
                {
                    ss << "\n" << indent;
                    accum = 0;
                }
            }
        }
        ss << "\n" << indent;
        for (auto& def : defs)
        {
            if (def.type == ArgType::unnamed)
            {
                const std::pair<std::string, std::string>& brace = braceType(def.requirement);
                std::string str                                  = brace.first + def.name + brace.second;
                ss << str << " ";
            }
        }
        ss << "\n";
        return ss.str();
    }
} // namespace codeUtils
