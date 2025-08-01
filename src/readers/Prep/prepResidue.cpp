#include "include/readers/Prep/prepResidue.hpp"

#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/readers/Prep/prepAtom.hpp"
#include "include/util/casting.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>

namespace gmml
{
    namespace prep
    {

        //////////////////////////////////////////////////////////
        //                       Constructor                    //
        //////////////////////////////////////////////////////////
        void initializePrepResidue(
            Residue* residue, PrepResidueProperties& properties, std::istream& in_file, std::string& line)
        {
            std::string dummy_atom_type;
            std::istringstream ss;
            properties.title = line; /// Set title of the residue
            getline(in_file, line);
            getline(in_file, line);
            // I guess you can extract with spaces as the default delimiter with >> from a stringstream
            ss.str(line);
            std::string name;
            ss >> name;
            residue->setName(name);
            properties.coordinateType = extractResidueCoordinateType(ss);
            properties.outputFormat = extractResidueOutputFormat(ss);
            ss.clear();
            getline(in_file, line);
            properties.geometryType = extractResidueGeometryType(ss);
            properties.dummyAtomOmission = extractResidueDummyAtomOmission(ss);
            ss >> dummy_atom_type;
            properties.dummyAtomType = dummy_atom_type;
            properties.dummyAtomPosition = extractResidueDummyAtomPosition(ss);
            getline(in_file, line);
            properties.charge = util::from_string<double>(line);
            /// Process atoms of the residue
            while (getline(in_file, line) && !util::Trim(line).empty())
            {
                Atom* atom = residue->addAtom(std::make_unique<Atom>());
                PrepAtomProperties atomProperties;
                initializePrepAtom(atom, atomProperties, line);
                properties.atomProperties.push_back(atomProperties);
            }
            /// Process the extra sections: IMPROPER, LOOP, DONE
            /// Skip blank lines until to reach to a known section title
            while (getline(in_file, line))
            { /// Does a corresponding action based on the section title
                switch (extractSectionType(line))
                {
                    case kSectionLoop:
                        properties.loops = extractLoops(in_file);
                        break;
                    case kSectionImproper:
                        properties.improperDihedrals = extractImproperDihedral(in_file);
                        break;
                    case kSectionDone:
                        return;
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
        }

        std::vector<std::string> getAtomNames(Residue* residue)
        {
            std::vector<std::string> foundAtoms;
            for (auto& prepAtom : residue->getAtoms())
            {
                foundAtoms.push_back(prepAtom->getName());
            }
            return foundAtoms;
        }

        std::vector<std::string> getHeavyAtomNames(Residue* residue)
        {
            std::vector<std::string> foundAtoms;
            for (auto& prepAtom : residue->getAtoms())
            {
                if (prepAtom->getName().at(0) != 'H' && prepAtom->getName() != "DUMM")
                {
                    foundAtoms.push_back(prepAtom->getName());
                }
            }
            return foundAtoms;
        }

        // Read about Amber Prep format. kTopTypeE is an "E" atom. The type determines connectivity and the order the
        // atoms are listed matters for what gets connected to what. This function goes through each atom and bonds it
        // to the correct atom according to the prep file tree structure. Loops are explicitly defined in their own
        // section. In Determine3dStructre if I travel up the first (because loops you can have two) incoming edge of
        // each atom I'll be traversing the atoms that define the bond, angle and dihedral so I can get a 3D coordinate
        // for each atom in order, starting with the DUMM atoms for the first real atom.
        void setConnectivities(Residue* residue, const PrepResidueProperties& properties)
        {
            std::vector<Atom*> atoms = residue->getAtoms();
            std::vector<uint> visitCount(atoms.size(), 0);
            std::vector<size_t> connectionPointStack {0};
            for (size_t n = 0; n < atoms.size(); n++)
            {
                size_t stackBack = connectionPointStack.back();
                addBond(atoms[stackBack], atoms[n]);
                visitCount[stackBack]++;
                if (visitCount[stackBack] >= properties.atomProperties[stackBack].topologicalType)
                {
                    connectionPointStack.pop_back();
                }
                if (properties.atomProperties[n].topologicalType > kTopTypeE)
                {
                    connectionPointStack.push_back(n);
                }
            }
            // Now bond any atoms defined in loops.
            std::vector<std::string> atomNames = gmml::atomNames(atoms);
            for (auto& loop : properties.loops)
            {
                size_t firstAtom = util::indexOf(atomNames, loop.first);
                size_t secondAtom = util::indexOf(atomNames, loop.second);
                if (firstAtom < atoms.size() && secondAtom < atoms.size())
                {
                    addBond(atoms[firstAtom], atoms[secondAtom]);
                }
            }
        }

        void generate3dStructure(Residue* residue, PrepResidueProperties& properties)
        { // Travel up the first incoming edge of each atom to traverse the atoms that define bond, angle and dihedral.
            // Each atoms position is defined by these connecting atoms. I set the first dummy atom position to 0,0,0.
            std::vector<Atom*> atoms = residue->getAtoms();
            if ((atoms.size() > 3) && (atoms.at(2)->getName() == "DUMM"))
            {
                // Set dummy atoms
                atoms[0]->setCoordinate(Coordinate(0, 0, 0));
                atoms[1]->setCoordinate(Coordinate(0.5, 0, 0));
                atoms[2]->setCoordinate(Coordinate(-0.75, 0.35, 0));
                // Use dummies as start for creating the other atoms.
                //		std::cout << "it1 is now pointing at atom: " << (*it1)->getName() << "\n";
                for (size_t n = 3; n < atoms.size(); n++)
                {
                    determine3dCoordinate(atoms[n], properties.atomProperties[n]);
                }
            }
            else // If we ever encounter the Kxyz type of prep file, write the code to handle it here.
            {
                std::string message = "Did not find dummy atoms in prep entry for " + residue->getName();
                util::log(__LINE__, __FILE__, util::ERR, message);
                throw std::runtime_error(message);
            }
            deleteDummyAtoms(residue, properties);
        }

        void deleteDummyAtoms(Residue* residue, PrepResidueProperties& properties)
        {
            std::vector<Atom*> atoms = residue->getAtoms();
            std::vector<std::string> names = atomNames(atoms);
            for (size_t n : util::reverse(util::indicesOfElement(names, std::string("DUMM"))))
            {
                residue->deleteAtom(atoms[n]);
                properties.atomProperties.erase(properties.atomProperties.begin() + n);
            }
        }

        /// Return coordinate type from a stream line which is the 2nd column of the 3rd line in each residue section
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

        /// Return geometry type from a stream line which is the first column of the 4th line in each residue section
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

        /// Parse the improper dihedral section of each residue section and return a std::vector of improper dihedrals
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
                if (atom_names[0].empty() || atom_names[1].empty() || atom_names[2].empty() || atom_names[3].empty())
                {
                    std::string message = "Error parsing DIHEDRAL section of prep file. Unexpected line >>>" + line +
                                          "<<<\nNote: DIHEDRAL sections must be bounded by empty lines!";
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

        //////////////////////////////////////////////////////////
        //                     FUNCTIONS                        //
        //////////////////////////////////////////////////////////
        double calculatePrepResidueCharge(Residue* residue)
        {
            double residue_charge = 0.0;
            for (auto& atom : residue->getAtoms())
            {
                residue_charge += atom->getCharge();
            }
            return residue_charge;
        }

        //////////////////////////////////////////////////////////
        //                     DISPLAY FUNCTIONS                //
        //////////////////////////////////////////////////////////
        std::string toString(Residue* residue, const PrepResidueProperties& properties)
        {
            std::stringstream out;
            out << "Title: " << properties.title << std::endl;
            out << std::setw(10) << "ResName" << std::setw(10) << "CrdType" << std::setw(10) << "Output"
                << std::setw(10) << "GeoType" << std::setw(15) << "DummyOmission" << std::setw(12) << "DummyType"
                << std::setw(12) << "DummyPos" << std::setw(10) << "Charge" << std::endl;
            out << std::setw(10) << residue->getName();

            if (properties.coordinateType == kINT)
            {
                out << std::setw(10) << "INT";
            }
            else if (properties.coordinateType == kXYZ)
            {
                out << std::setw(10) << "XYZ";
            }
            else
            {
                out << std::setw(10) << "--";
            }

            if (properties.outputFormat == kBinary)
            {
                out << std::setw(10) << "Binary";
            }
            else if (properties.outputFormat == kFormatted)
            {
                out << std::setw(10) << "NBinary";
            }
            else
            {
                out << std::setw(10) << "--";
            }

            if (properties.geometryType == kGeometryCorrect)
            {
                out << std::setw(10) << "Correct";
            }
            else if (properties.geometryType == kGeometryChange)
            {
                out << std::setw(10) << "Change";
            }
            else
            {
                out << std::setw(10) << "--";
            }

            if (properties.dummyAtomOmission == kOmit)
            {
                out << std::setw(15) << "YES";
            }
            else if (properties.dummyAtomOmission == kNomit)
            {
                out << std::setw(15) << "NO";
            }
            else
            {
                out << std::setw(15) << "--";
            }

            out << std::setw(12) << properties.dummyAtomType;

            if (properties.dummyAtomPosition == kPositionAll)
            {
                out << std::setw(12) << "ALL";
            }
            else if (properties.dummyAtomPosition == kPositionBeg)
            {
                out << std::setw(12) << "BEG";
            }
            else
            {
                out << std::setw(12) << "--";
            }

            out << std::setw(10) << properties.charge << std::endl << std::endl;

            out << std::setw(3) << "#" << std::setw(6) << "Name" << std::setw(6) << "Type" << std::setw(3) << "TT"
                << std::setw(4) << "B#" << std::setw(4) << "A#" << std::setw(4) << "D#" << std::setw(10) << "Bond"
                << std::setw(10) << "Angle" << std::setw(10) << "Dihedral" << std::setw(10) << "Charge" << std::setw(10)
                << "Bonded" << std::endl;

            out << std::endl << "Improper dihedrals" << std::endl;
            for (auto& improperDihedral : properties.improperDihedrals)
            {
                for (auto& dihedralComponent : improperDihedral)
                {
                    out << std::setw(6) << dihedralComponent;
                }
                out << std::endl;
            }

            out << std::endl << "Loops" << std::endl;

            out << std::endl;
            return out.str();
        }

        void write(Residue* residue, const PrepResidueProperties& properties, std::ostream& stream)
        {
            stream << properties.title << std::endl
                   << std::endl
                   << std::left << std::setw(4) << residue->getName() << " " << std::right << std::setw(3)
                   << coordinateTypeNames[properties.coordinateType] << " " << std::setw(1) << properties.outputFormat
                   << std::endl
                   << geometryTypesNames[properties.geometryType] << " "
                   << dummyAtomOmissionNames[properties.dummyAtomOmission] << " " << properties.dummyAtomType << " "
                   << dummyAtomPositionNames[properties.dummyAtomPosition] << std::endl
                   << std::right << std::setw(8) << std::fixed << std::setprecision(3) << properties.charge
                   << std::endl;
            std::vector<Atom*> atoms = residue->getAtoms();
            for (size_t n = 0; n < atoms.size(); n++)
            {
                write(atoms[n], properties.atomProperties[n], stream);
            }
            stream << std::endl;
            if (properties.improperDihedrals.size() > 0)
            {
                stream << "IMPROPER" << std::endl;
                for (auto& dihedral : properties.improperDihedrals)
                {
                    stream << std::left << std::setw(4) << dihedral.at(0) << std::left << std::setw(4) << dihedral.at(1)
                           << std::left << std::setw(4) << dihedral.at(2) << std::left << std::setw(4) << dihedral.at(3)
                           << std::endl;
                }
                stream << std::endl;
            }
            if (properties.loops.size() > 0)
            {
                stream << "LOOP" << std::endl;
                for (auto& loop : properties.loops)
                {
                    stream << loop.first << " " << loop.second << std::endl;
                }
            }
            stream << std::endl << "DONE" << std::endl;
        }
    } // namespace prep
} // namespace gmml
