#ifndef INCLUDE_UTIL_STRUCTUREDFILES_HPP
#define INCLUDE_UTIL_STRUCTUREDFILES_HPP

#include <ostream>
#include <string>
#include <variant>
#include <vector>

namespace gmml
{
    namespace util
    {
        struct TextHeader
        {
            uint size;
            std::string text;
        };

        struct TextParagraph
        {
            std::vector<std::string> lines;
        };

        struct TextTable
        {
            std::vector<std::string> header;
            std::vector<std::vector<std::string>> rows;
        };

        typedef std::variant<TextHeader, TextParagraph, TextTable> TextVariant;

        void toTxt(std::ostream& stream, const std::vector<TextVariant>& contents);
        void toHtml(std::ostream& stream, const std::vector<TextVariant>& contents);
        void toCsv(std::ostream& stream, const std::string& delimiter, const TextTable& table);
    } // namespace util
} // namespace gmml
#endif
