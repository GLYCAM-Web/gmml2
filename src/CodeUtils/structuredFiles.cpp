#include "includes/CodeUtils/structuredFiles.hpp"

#include <string>
#include <vector>
#include <variant>
#include <ostream>

namespace codeUtils
{
    void toTxt(std::ostream& stream, const std::vector<TextVariant>& contents)
    {
        auto printHeader = [&stream](const TextHeader& a)
        {
            uint hashCount = std::max(uint(0), 7 - a.size);
            std::string hashes(hashCount, '#');
            stream << hashes << " " << a.text << " " << hashes << "\n\n";
        };

        auto printParagraph = [&stream](const TextParagraph& a)
        {
            for (auto& line : a.lines)
            {
                stream << line << "\n";
            }
            stream << "\n";
        };

        auto printTable = [&stream](const TextTable& table)
        {
            size_t columns = table.header.size();
            std::vector<size_t> columnSize(columns, 0);
            for (size_t n = 0; n < columns; n++)
            {
                size_t size = table.header[n].size();
                for (auto& row : table.rows)
                {
                    size = std::max(size, row[n].size());
                }
                columnSize[n] = size;
            }
            auto printRow = [&](const std::string& delimiter, const std::vector<std::string>& row)
            {
                for (size_t n = 0; n < columns; n++)
                {
                    const std::string& str = row[n];
                    stream << str << std::string(columnSize[n] - str.size(), ' ');
                    if (n < columns - 1)
                    {
                        stream << delimiter;
                    }
                }
                stream << "\n";
            };
            printRow(" | ", table.header);
            for (auto& row : table.rows)
            {
                printRow("   ", row);
            }
            stream << "\n";
        };

        for (auto& a : contents)
        {
            if (std::holds_alternative<TextHeader>(a))
            {
                printHeader(std::get<TextHeader>(a));
            }
            else if (std::holds_alternative<TextParagraph>(a))
            {
                printParagraph(std::get<TextParagraph>(a));
            }
            else if (std::holds_alternative<TextTable>(a))
            {
                printTable(std::get<TextTable>(a));
            }
            else
            {
                throw std::runtime_error("unhandled variant in toTxt");
            }
        }
    }

    void toHtml(std::ostream& stream, const std::vector<TextVariant>& contents)
    {
        auto printHeader = [&stream](const TextHeader& a)
        {
            std::string h = "h" + std::to_string(a.size);
            stream << "<" << h << ">" << a.text << "</" << h << ">\n";
        };

        auto printParagraph = [&stream](const TextParagraph& a)
        {
            for (auto& line : a.lines)
            {
                stream << "<p>" << line << "</p>\n";
            }
        };

        auto printTable = [&stream](const codeUtils::TextTable& table)
        {
            auto printRow = [&](const std::vector<std::string>& row)
            {
                stream << "<tr>\n";
                for (auto& str : row)
                {
                    stream << "<th>" << str << "</th>\n";
                }
                stream << "</tr>\n";
            };
            stream << "<table>\n";
            printRow(table.header);
            for (auto& row : table.rows)
            {
                printRow(row);
            }
            stream << "</table>\n";
        };

        stream << "<!DOCTYPE html>\n"
               << "<html>\n"
               << "<head>\n"
               << "<style>\n"
               << "table {\n"
               << "font-family: arial, sans-serif;\n"
               << "border-collapse: collapse;\n"
               << "}\n"
               << "td, th {\n"
               << "border: 1px solid #dddddd;\n"
               << "text-align: left;\n"
               << "padding: 8px;\n"
               << "}\n"
               << "tr:nth-child(even) {\n"
               << "background-color: #dddddd;\n"
               << "}\n"
               << "</style>\n"
               << "</head>\n"
               << "<body>\n";

        for (auto& a : contents)
        {
            if (std::holds_alternative<TextHeader>(a))
            {
                printHeader(std::get<TextHeader>(a));
            }
            else if (std::holds_alternative<TextParagraph>(a))
            {
                printParagraph(std::get<TextParagraph>(a));
            }
            else if (std::holds_alternative<TextTable>(a))
            {
                printTable(std::get<TextTable>(a));
            }
            else
            {
                throw std::runtime_error("unhandled variant in toHtml");
            }
        }

        stream << "</body>\n"
               << "</html>\n";
    }

    void toCsv(std::ostream& stream, const std::string& delimiter, const TextTable& table)
    {
        auto printRow = [&](const std::string& linebreak, const std::vector<std::string>& row)
        {
            for (size_t n = 0; n < row.size(); n++)
            {
                stream << row[n] << (n < row.size() - 1 ? delimiter : linebreak);
            }
        };
        printRow("\n", table.header);
        for (size_t n = 0; n < table.rows.size(); n++)
        {
            std::string linebreak = (n == table.rows.size() - 1) ? "" : "\n";
            printRow(linebreak, table.rows[n]);
        }
    }
} // namespace codeUtils
