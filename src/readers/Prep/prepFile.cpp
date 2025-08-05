#include "include/readers/Prep/prepFile.hpp"

#include "include/readers/Prep/prepAtom.hpp"
#include "include/readers/Prep/prepResidue.hpp"
#include "include/util/casting.hpp"
#include "include/util/containers.hpp"
#include "include/util/files.hpp"
#include "include/util/filesystem.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <ios>

namespace gmml
{
    namespace prep
    {
        PrepData readPrepFile(const std::string& prep_file)
        {
            PrepData data;
            util::ensureFileExists(prep_file);
            std::ifstream in_file(prep_file.c_str());
            if (in_file.is_open())
            {
                readAllResidues(data, in_file);
                in_file.close();
            }
            else
            {
                throw std::runtime_error("Prep file exists but couldn't be opened.");
            }
            setAtomConnectivities(data);
            generate3dStructures(data);
            return data;
        }

        PrepData readPrepFile(const std::string& prep_file, const std::vector<std::string> queryNames)
        {
            PrepData data;
            util::ensureFileExists(prep_file);
            std::ifstream in_file(prep_file.c_str());
            if (in_file.is_open())
            {
                readQueryResidues(data, in_file, queryNames);
                in_file.close();
            }
            else
            {
                throw std::runtime_error("Prep file exists but couldn't be opened.");
            }
            // Note that I'm assuming that given a list of query names, the user wants
            // to use all these residues, thus this isn't wasteful:
            util::log(__LINE__, __FILE__, util::INF, "Finished reading prep file. Now setting atomic connectivities.");
            setAtomConnectivities(data);
            util::log(
                __LINE__, __FILE__, util::INF, "Finished setting atomic connectivities. Now generating 3D structures.");
            generate3dStructures(data);
            util::log(__LINE__, __FILE__, util::INF, "Finished, returning from PreFile constructor..");
            return data;
        }

        //////////////////////////////////////////////////////////
        //                         MUTATORS                     //
        //////////////////////////////////////////////////////////
        void setAtomConnectivities(PrepData& data)
        {
            for (size_t n = 0; n < data.residueCount; n++)
            {
                setConnectivities(data, n);
            }
        }

        void generate3dStructures(PrepData& data)
        {
            for (size_t n = 0; n < data.residueCount; n++)
            {
                generate3dStructure(data, n);
            }
        }

        //////////////////////////////////////////////////////////
        //                         FUNCTIONS                    //
        //////////////////////////////////////////////////////////
        void readAllResidues(PrepData& data, std::istream& in_file)
        {
            std::string line = "";
            getline(in_file, line);
            getline(in_file, line); // first two lines are always blank apparently. smh.
            getline(in_file, line); // This should be first line of residue entry. Title
            while (util::Trim(line).find("STOP") == std::string::npos) /// End of file
            {
                initializePrepResidue(data, in_file, line);
                getline(in_file, line); // This should be first line of next residue entry or STOP.
            }
        }

        // Reads each line of the file. If it finds one of the query residues it reads it in. Won't read in query
        // repeats twice.
        void readQueryResidues(PrepData& data, std::istream& in_file, const std::vector<std::string>& queryNames)
        {
            std::string line = "";
            getline(in_file, line);
            getline(in_file, line); // first two lines of the file are always blank apparently. smh.
            getline(in_file, line); // This should be first line of residue entry. Title.
            while (util::Trim(line).find("STOP") == std::string::npos) // While not at end of file
            {
                std::streampos firstResidueLinePosition =
                    in_file.tellg(); // save correct position to start reading residue
                std::string savedTitle = line;
                // Need to move to line with residue name on it.
                getline(in_file, line);                                            // blank line
                getline(in_file, line);                                            // residue name appears here
                std::vector<std::string> residueNameLine = util::split(line, ' '); // front() string will be name
                if (util::contains(queryNames, residueNameLine.front()))
                {
                    int numberOfTimesToReadInResidue =
                        std::count(queryNames.begin(), queryNames.end(), residueNameLine.front());
                    while (numberOfTimesToReadInResidue > 0)
                    {
                        in_file.seekg(firstResidueLinePosition); // go back here so the residue constructor works
                        line = savedTitle;
                        initializePrepResidue(data, in_file, line);
                        --numberOfTimesToReadInResidue;
                    }
                }
                else
                { // need to flush the lines until we find a residue we want.
                    while (util::Trim(line).find("DONE") == std::string::npos)
                    {
                        getline(in_file, line);
                    }
                }
                getline(in_file, line); // This should be first line of next residue entry or STOP.
            }
        }

        void write(const PrepData& data, const std::string& prep_file)
        {
            try
            {
                util::writeToFile(prep_file, [&](std::ostream& stream) { write(data, stream); });
            }
            catch (...)
            {
                throw std::runtime_error("PrepFile could not be created for writing");
            }
        }

        void write(const PrepData& data, std::ostream& stream)
        {
            stream << "\n"
                   << "\n";
            for (size_t n = 0; n < data.residueCount; n++)
            {
                writeResidue(data, n, stream);
            }
            stream << "STOP\n";
        }

        std::string print(const PrepData& data)
        {
            std::string out;
            for (size_t n = 0; n < data.residueCount; n++)
            {
                out += "**********************************************************************************\n";
                out += residueToString(data, n);
            }
            return out;
        }
    } // namespace prep
} // namespace gmml
