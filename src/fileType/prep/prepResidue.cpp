#include "include/fileType/prep/prepResidue.hpp"

#include "include/fileType/prep/prepAtom.hpp"
#include "include/fileType/prep/prepFunctions.hpp"
#include "include/graph/graphManipulation.hpp"
#include "include/util/casting.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <fstream>
#include <iomanip>
#include <ios>
#include <sstream>

namespace gmml
{
    namespace prep
    {
        namespace
        {
            /// Return coordinate type from a stream line which is the 2nd column of the 3rd line in each residue
            /// section
            CoordinateType extractResidueCoordinateType(std::istream& ss)
            {
                std::string s;
                ss >> s;
                return (s == "XYZ") ? kXYZ : kINT;
            }

            /// Return output format from a stream line which is the 3rd column of the 3rd line in each residue section
            OutputFormat extractResidueOutputFormat(std::istream& ss)
            {
                int val;
                ss >> val;
                return (val == 1) ? kBinary : kFormatted;
            }

            /// Return geometry type from a stream line which is the first column of the 4th line in each residue
            /// section
            GeometryType extractResidueGeometryType(std::istream& ss)
            {
                std::string s;
                ss >> s;
                return (s == "CHANGE") ? kGeometryChange : kGeometryCorrect;
            }

            /// Return dummy atom omission from a stream line which is the 2nd column of the 4th line in each residue
            /// section
            DummyAtomOmission extractResidueDummyAtomOmission(std::istream& ss)
            {
                std::string s;
                ss >> s;
                return (s == "NOMIT") ? kNomit : kOmit;
            }

            /// Return dummy atom position from a stream line which is the 4th column of the 4th line in each residue
            /// section
            DummyAtomPosition extractResidueDummyAtomPosition(std::istream& ss)
            {
                std::string s;
                ss >> s;
                return (s == "ALL") ? kPositionAll : kPositionBeg;
            }

            /// Return a corresponding title from a stream line which may appear in each residue section
            SectionType extractSectionType(std::string& line)
            {
                if (line == "LOOP")
                {
                    return kSectionLoop;
                }
                else if (line == "IMPROPER")
                {
                    return kSectionImproper;
                }
                else if (line == "DONE")
                {
                    return kSectionDone;
                }
                return kSectionOther;
            }

            /// Parse the loop section of each residue section and return a loop map
            std::vector<std::pair<std::string, std::string>> extractLoops(std::istream& in_file)
            {
                std::vector<std::pair<std::string, std::string>> result;
                std::string line;
                std::stringstream ss;
                while (getline(in_file, line) &&
                       !util::Trim(line).empty()) /// Read file until blank line which determines the end of the section
                {
                    ss.clear();
                    ss.str(line);
                    std::string atom_names[2];
                    ss >> atom_names[0] >> atom_names[1];
                    if (atom_names[0].empty() || atom_names[1].empty())
                    {
                        std::string message = "Error parsing LOOP section of prep file. Unexpected line >>>" + line +
                                              "<<<\nNote: LOOP sections must be bounded by empty lines!";
                        util::log(__LINE__, __FILE__, util::ERR, message);
                        throw std::runtime_error(message);
                    }
                    result.push_back({atom_names[0], atom_names[1]});
                }
                return result;
            }

            /// Parse the improper dihedral section of each residue section and return a std::vector of improper
            /// dihedrals
            DihedralVector extractImproperDihedral(std::istream& in_file)
            {
                std::string line;
                std::stringstream ss;
                DihedralVector dihedrals;
                while (getline(in_file, line) &&
                       !util::Trim(line).empty()) /// Read file until blank line which determines the end of the section
                {
                    std::string atom_names[4];
                    Dihedral dihedral;
                    ss.clear();
                    ss.str(line);
                    /// Extract improper atom types involving in a dihedral from each line
                    ss >> atom_names[0] >> atom_names[1] >> atom_names[2] >> atom_names[3];
                    if (atom_names[0].empty() || atom_names[1].empty() || atom_names[2].empty() ||
                        atom_names[3].empty())
                    {
                        std::string message = "Error parsing DIHEDRAL section of prep file. Unexpected line >>>" +
                                              line + "<<<\nNote: DIHEDRAL sections must be bounded by empty lines!";
                        util::log(__LINE__, __FILE__, util::ERR, message);
                        throw std::runtime_error(message);
                    }
                    for (int i = 0; i < 4; i++)
                    {
                        dihedral.push_back(atom_names[i]);
                    }
                    dihedrals.push_back(dihedral); /// Create a new dihedral into the std::vector of dihedrals
                }
                return dihedrals;
            }
        } // namespace

        void initializePrepResidue(PrepData& data, std::istream& in_file, std::string& line)
        {
            std::string dummy_atom_type;
            std::istringstream ss;
            size_t residueId = data.residueCount;
            data.residueCount++;
            data.residues.title.push_back(line);
            getline(in_file, line);
            getline(in_file, line);
            ss.str(line);
            std::string name;
            ss >> name;
            data.residues.name.push_back(name);
            data.residues.coordinateType.push_back(extractResidueCoordinateType(ss));
            data.residues.outputFormat.push_back(extractResidueOutputFormat(ss));
            ss.clear();
            getline(in_file, line);
            data.residues.geometryType.push_back(extractResidueGeometryType(ss));
            data.residues.dummyAtomOmission.push_back(extractResidueDummyAtomOmission(ss));
            ss >> dummy_atom_type;
            data.residues.dummyAtomType.push_back(dummy_atom_type);
            data.residues.dummyAtomPosition.push_back(extractResidueDummyAtomPosition(ss));
            getline(in_file, line);
            data.residues.charge.push_back(util::from_string<double>(line));

            while (getline(in_file, line) && !util::Trim(line).empty())
            {
                initializePrepAtom(data, residueId, line);
            }
            /// Process the extra sections: IMPROPER, LOOP, DONE
            /// Skip blank lines until to reach to a known section title

            DihedralVector improperDihedrals;
            std::vector<std::pair<std::string, std::string>> loops;
            bool keepGoing = true;
            while (keepGoing && getline(in_file, line))
            { /// Does a corresponding action based on the section title
                switch (extractSectionType(line))
                {
                    case kSectionLoop:
                        loops = extractLoops(in_file);
                        break;
                    case kSectionImproper:
                        improperDihedrals = extractImproperDihedral(in_file);
                        break;
                    case kSectionDone:
                        keepGoing = false;
                        break;
                    case kSectionBlank:
                        break;
                    case kSectionOther:
                        util::log(
                            __LINE__,
                            __FILE__,
                            util::WAR,
                            "Unrecognized section in prep file inside these brackets >>>" + line + "<<<");
                        break;
                    default:
                        util::log(
                            __LINE__,
                            __FILE__,
                            util::WAR,
                            "Did not figure out a prep section type for line inside these brackets >>>" + line + "<<<");
                }
            }
            data.residues.loops.push_back(loops);
            data.residues.improperDihedrals.push_back(improperDihedrals);
        }

        std::vector<std::string> getAtomNames(PrepData& data, size_t residueId)
        {
            return util::indicesToValues(data.atoms.name, residueAtoms(data, residueId));
        }

        std::vector<std::string> getHeavyAtomNames(PrepData& data, size_t residueId)
        {
            std::function<bool(const std::string&)> isHeavy = [](const std::string& str)
            { return !(str.empty() || str[0] == 'H' || str == "DUMM"); };

            return util::vectorFilter(isHeavy, getAtomNames(data, residueId));
        }

        // Read about Amber Prep format. kTopTypeE is an "E" atom. The type determines connectivity and the order the
        // atoms are listed matters for what gets connected to what. This function goes through each atom and bonds it
        // to the correct atom according to the prep file tree structure. Loops are explicitly defined in their own
        // section. In Determine3dStructre if I travel up the first (because loops you can have two) incoming edge of
        // each atom I'll be traversing the atoms that define the bond, angle and dihedral so I can get a 3D coordinate
        // for each atom in order, starting with the DUMM atoms for the first real atom.
        void setConnectivities(PrepData& data, size_t index)
        {
            std::vector<size_t> atomIds = residueAtoms(data, index);
            std::vector<uint> visitCount(data.atomCount, 0);
            std::vector<size_t> connectionPointStack {atomIds[0]};
            for (size_t n : atomIds)
            {
                size_t stackBack = connectionPointStack.back();
                if (stackBack != n)
                {
                    graph::addEdge(data.atomGraph, {stackBack, n});
                }
                visitCount[stackBack]++;
                if (visitCount[stackBack] >= data.atoms.topologicalType[stackBack])
                {
                    connectionPointStack.pop_back();
                }
                if (data.atoms.topologicalType[n] > kTopTypeE)
                {
                    connectionPointStack.push_back(n);
                }
            }
            // Now bond any atoms defined in loops.
            std::vector<std::string> atomNames = util::indicesToValues(data.atoms.name, atomIds);
            for (auto& loop : data.residues.loops[index])
            {
                size_t firstAtom = util::indexOf(atomNames, loop.first);
                size_t secondAtom = util::indexOf(atomNames, loop.second);
                if (firstAtom < atomIds.size() && secondAtom < atomIds.size())
                {
                    graph::addEdge(data.atomGraph, {atomIds[firstAtom], atomIds[secondAtom]});
                }
            }
        }

        void generate3dStructure(PrepData& data, size_t index)
        { // Travel up the first incoming edge of each atom to traverse the atoms that define bond, angle and dihedral.
            // Each atoms position is defined by these connecting atoms. I set the first dummy atom position to 0,0,0.
            std::vector<size_t> atomIds = residueAtoms(data, index);
            if ((atomIds.size() > 3) && (data.atoms.name[atomIds[2]] == "DUMM"))
            {
                // Set dummy atoms
                std::vector<Coordinate> coords {
                    {    0,    0, 0},
                    {  0.5,    0, 0},
                    {-0.75, 0.35, 0}
                };
                for (size_t k = 0; k < 3; k++)
                {
                    data.atoms.coordinate[atomIds[k]] = coords[k];
                }
                // Use dummies as start for creating the other atoms.
                for (size_t n = 3; n < atomIds.size(); n++)
                {
                    data.atoms.coordinate[atomIds[n]] = determine3dCoordinate(data, atomIds[n]);
                }
            }
            else // If we ever encounter the Kxyz type of prep file, write the code to handle it here.
            {
                std::string message = "Did not find dummy atoms in prep entry for " + data.residues.name[index];
                util::log(__LINE__, __FILE__, util::ERR, message);
                throw std::runtime_error(message);
            }
            deleteDummyAtoms(data, index);
        }

        void deleteDummyAtoms(PrepData& data, size_t index)
        {
            std::vector<size_t> atomIds = residueAtoms(data, index);
            std::vector<std::string> names = util::indicesToValues(data.atoms.name, atomIds);
            for (size_t k : util::indicesOfElement(names, std::string("DUMM")))
            {
                graph::removeNode(data.atomGraph, atomIds[k]);
            }
        }

        //////////////////////////////////////////////////////////
        //                     FUNCTIONS                        //
        //////////////////////////////////////////////////////////
        double calculatePrepResidueCharge(PrepData& data, size_t residueId)
        {
            return util::vectorSum(0.0, util::indicesToValues(data.atoms.charge, residueAtoms(data, residueId)));
        }

        //////////////////////////////////////////////////////////
        //                     DISPLAY FUNCTIONS                //
        //////////////////////////////////////////////////////////
        std::string residueToString(const PrepData& data, size_t index)
        {
            const ResidueData& residues = data.residues;
            std::stringstream out;
            out << "Title: " << residues.title[index] << std::endl
                << std::setw(10) << "ResName" << std::setw(10) << "CrdType" << std::setw(10) << "Output"
                << std::setw(10) << "GeoType" << std::setw(15) << "DummyOmission" << std::setw(12) << "DummyType"
                << std::setw(12) << "DummyPos" << std::setw(10) << "Charge" << std::endl
                << std::setw(10) << residues.name[index] << std::setw(10)
                << coordinateTypeNames[residues.coordinateType[index]] << std::setw(10)
                << outputFormatNames[residues.outputFormat[index]] << std::setw(10)
                << geometryTypesNames[residues.geometryType[index]] << std::setw(15)
                << dummyAtomOmissionYesNo[residues.dummyAtomOmission[index]] << std::setw(12)
                << residues.dummyAtomType[index] << std::setw(12)
                << dummyAtomPositionNames[residues.dummyAtomPosition[index]] << std::setw(10) << residues.charge[index]
                << std::endl
                << std::endl
                << std::setw(3) << "#" << std::setw(6) << "Name" << std::setw(6) << "Type" << std::setw(3) << "TT"
                << std::setw(4) << "B#" << std::setw(4) << "A#" << std::setw(4) << "D#" << std::setw(10) << "Bond"
                << std::setw(10) << "Angle" << std::setw(10) << "Dihedral" << std::setw(10) << "Charge" << std::setw(10)
                << "Bonded" << std::endl
                << std::endl
                << "Improper dihedrals" << std::endl;
            for (auto& improperDihedral : residues.improperDihedrals[index])
            {
                for (auto& dihedralComponent : improperDihedral)
                {
                    out << std::setw(6) << dihedralComponent;
                }
                out << std::endl;
            }

            out << std::endl << "Loops" << std::endl << std::endl;
            return out.str();
        }

        void writeResidue(const PrepData& data, size_t index, std::ostream& stream)
        {
            const ResidueData& residues = data.residues;
            stream << residues.title[index] << std::endl
                   << std::endl
                   << std::left << std::setw(4) << residues.name[index] << " " << std::right << std::setw(3)
                   << coordinateTypeNames[residues.coordinateType[index]] << " " << std::setw(1)
                   << residues.outputFormat[index] << std::endl
                   << geometryTypesNames[residues.geometryType[index]] << " "
                   << dummyAtomOmissionNames[residues.dummyAtomOmission[index]] << " " << residues.dummyAtomType[index]
                   << " " << dummyAtomPositionNames[residues.dummyAtomPosition[index]] << std::endl
                   << std::right << std::setw(8) << std::fixed << std::setprecision(3) << residues.charge[index]
                   << std::endl;
            std::vector<size_t> atomIds = residueAtoms(data, index);
            size_t atomCount = atomIds.size();
            for (size_t n = 0; n < atomCount; n++)
            {
                writeAtom(data, atomIds[n], stream);
            }
            stream << std::endl;
            if (residues.improperDihedrals[index].size() > 0)
            {
                stream << "IMPROPER" << std::endl;
                for (auto& dihedral : residues.improperDihedrals[index])
                {
                    stream << std::left << std::setw(4) << dihedral.at(0) << std::left << std::setw(4) << dihedral.at(1)
                           << std::left << std::setw(4) << dihedral.at(2) << std::left << std::setw(4) << dihedral.at(3)
                           << std::endl;
                }
                stream << std::endl;
            }
            if (residues.loops[index].size() > 0)
            {
                stream << "LOOP" << std::endl;
                for (auto& loop : residues.loops[index])
                {
                    stream << loop.first << " " << loop.second << std::endl;
                }
            }
            stream << std::endl << "DONE" << std::endl;
        }
    } // namespace prep
} // namespace gmml
