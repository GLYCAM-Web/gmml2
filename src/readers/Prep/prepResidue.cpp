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
        PrepResidue::PrepResidue(std::istream& in_file, std::string& line)
        {
            std::string dummy_atom_type;
            std::istringstream ss;
            properties.title = line; /// Set title of the residue
            getline(in_file, line);
            getline(in_file, line);
            // I guess you can extract with spaces as the default delimiter with >> from a stringstream
            ss.str(line);
            this->ExtractResidueName(ss);
            this->ExtractResidueCoordinateType(ss);
            this->ExtractResidueOutputFormat(ss);
            ss.clear();
            getline(in_file, line);
            this->ExtractResidueGeometryType(ss);
            this->ExtractResidueDummyAtomOmission(ss);
            ss >> dummy_atom_type;
            properties.dummyAtomType = dummy_atom_type;
            this->ExtractResidueDummyAtomPosition(ss);
            getline(in_file, line);
            properties.charge = util::from_string<double>(line);
            /// Process atoms of the residue
            while (getline(in_file, line) && !util::Trim(line).empty())
            {
                this->addAtom(std::make_unique<PrepAtom>(line));
            }
            /// Process the extra sections: IMPROPER, LOOP, DONE
            /// Skip blank lines until to reach to a known section title
            while (getline(in_file, line))
            { /// Does a corresponding action based on the section title
                switch (ExtractSectionType(line))
                {
                    case kSectionLoop:
                        this->ExtractLoops(in_file);
                        break;
                    case kSectionImproper:
                        this->ExtractImproperDihedral(in_file);
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

        std::string PrepResidue::GetStringFormatOfCoordinateType(CoordinateType coordinate_type) const
        {
            switch (coordinate_type)
            {
                case kINT:
                    return "INT";
                case kXYZ:
                    return "XYZ";
                default:
                    return "";
            }
        }

        std::string PrepResidue::GetStringFormatOfGeometryType(GeometryType geometry_type) const
        {
            switch (geometry_type)
            {
                case kGeometryCorrect:
                    return "CORRECT";
                case kGeometryChange:
                    return "CHANGE";
                default:
                    return "";
            }
        }

        std::string PrepResidue::GetStringFormatOfDummyAtomPosition(DummyAtomPosition dummy_atom_position) const
        {
            switch (dummy_atom_position)
            {
                case kPositionAll:
                    return "ALL";
                case kPositionBeg:
                    return "BEG";
                default:
                    return "";
            }
        }

        std::string PrepResidue::GetStringFormatOfDummyAtomOmission(DummyAtomOmission dummy_atom_omission) const
        {
            switch (dummy_atom_omission)
            {
                case kOmit:
                    return "OMIT";
                case kNomit:
                    return "NOMIT";
                default:
                    return "";
            }
        }

        std::vector<std::string> PrepResidue::GetAtomNames() const
        {
            std::vector<std::string> foundAtoms;
            for (auto& prepAtom : this->getAtoms())
            {
                foundAtoms.push_back(prepAtom->getName());
            }
            return foundAtoms;
        }

        std::vector<std::string> PrepResidue::GetHeavyAtomNames() const
        {
            std::vector<std::string> foundAtoms;
            for (auto& prepAtom : this->getAtoms())
            {
                if (prepAtom->getName().at(0) != 'H' && prepAtom->getName() != "DUMM")
                {
                    foundAtoms.push_back(prepAtom->getName());
                }
            }
            return foundAtoms;
        }

        //////////////////////////////////////////////////////////
        //                         FUNCTIONS                    //
        //////////////////////////////////////////////////////////
        // Read about Amber Prep format. kTopTypeE is an "E" atom. The type determines connectivity and the order the
        // atoms are listed matters for what gets connected to what. This function goes through each atom and bonds it
        // to the correct atom according to the prep file tree structure. Loops are explicitly defined in their own
        // section. In Determine3dStructre if I travel up the first (because loops you can have two) incoming edge of
        // each atom I'll be traversing the atoms that define the bond, angle and dihedral so I can get a 3D coordinate
        // for each atom in order, starting with the DUMM atoms for the first real atom.
        void PrepResidue::SetConnectivities()
        {
            // std::cout << "Attempting to set connectivities in prepResidue: " << this->getName() << std::endl;
            // std::cout << "Number of atoms: " << this->getAtoms().size() << std::endl;
            // std::cout << "First atom is " << this->getAtoms().front()->getName() << std::endl;
            std::vector<PrepAtom*> connectionPointStack = {util::erratic_cast<PrepAtom*>(this->getAtoms().front())};
            // while(currentAtom != this->getAtoms().end())
            for (auto& currentAtom : this->getAtoms())
            {
                PrepAtom* currentAtomAsPrepType = util::erratic_cast<PrepAtom*>(currentAtom);
                addBond(connectionPointStack.back(), currentAtomAsPrepType);
                // std::cout << "Bonded " << connectionPointStack.back()->getName() << " to " << currentAtom->getName()
                // << std::endl;;
                connectionPointStack.back()->properties.visitCount++;
                if (connectionPointStack.back()->properties.visitCount >=
                    connectionPointStack.back()->properties.topologicalType)
                {
                    connectionPointStack.pop_back();
                }
                if (currentAtomAsPrepType->properties.topologicalType > kTopTypeE)
                {
                    connectionPointStack.push_back(currentAtomAsPrepType);
                }
                ++currentAtom;
            }
            // Now bond any atoms defined in loops.
            std::vector<Atom*> atoms = this->getAtoms();
            std::vector<std::string> atomNames = gmml::atomNames(atoms);
            for (auto& loop : this->properties.loops)
            {
                size_t firstAtom = util::indexOf(atomNames, loop.first);
                size_t secondAtom = util::indexOf(atomNames, loop.second);
                if (firstAtom < atoms.size() && secondAtom < atoms.size())
                {
                    addBond(atoms[firstAtom], atoms[secondAtom]);
                }
            }
        }

        void PrepResidue::Generate3dStructure()
        { // Travel up the first incoming edge of each atom to traverse the atoms that define bond, angle and dihedral.
            // Each atoms position is defined by these connecting atoms. I set the first dummy atom position to 0,0,0.
            if ((this->getAtoms().size() > 3) && (this->getAtoms().at(2)->getName() == "DUMM"))
            {
                // Set dummy atoms
                this->getAtoms().at(0)->setCoordinate(Coordinate(0, 0, 0));
                this->getAtoms().at(1)->setCoordinate(Coordinate(0.5, 0, 0));
                this->getAtoms().at(2)->setCoordinate(Coordinate(-0.75, 0.35, 0));
                // Use dummies as start for creating the other atoms.
                std::vector<Atom*> atomsInResidue = this->getAtoms();
                std::vector<Atom*>::iterator it1 = atomsInResidue.begin();
                std::advance(it1, 3);
                //		std::cout << "it1 is now pointing at atom: " << (*it1)->getName() << "\n";
                while (it1 != atomsInResidue.end())
                {
                    PrepAtom* it1AsPrepAtom = util::erratic_cast<PrepAtom*>(*it1);
                    it1AsPrepAtom->Determine3dCoordinate();
                    ++it1;
                }
            }
            else // If we ever encounter the Kxyz type of prep file, write the code to handle it here.
            {
                std::string message = "Did not find dummy atoms in prep entry for " + this->getName();
                util::log(__LINE__, __FILE__, util::ERR, message);
                throw std::runtime_error(message);
            }
            this->DeleteDummyAtoms();
        }

        void PrepResidue::DeleteDummyAtoms()
        {
            std::vector<Atom*> atoms = this->getAtoms();
            std::vector<std::string> names = atomNames(atoms);
            for (size_t n : util::reverse(util::indicesOfElement(names, std::string("DUMM"))))
            {
                this->deleteAtom(atoms[n]);
            }
        }

        /// Return residue name from a stream line which is the first column of the 3rd line in each residue section
        void PrepResidue::ExtractResidueName(std::istream& ss)
        {
            std::string name;
            ss >> name;
            this->setName(name);
        }

        /// Return coordinate type from a stream line which is the 2nd column of the 3rd line in each residue section
        void PrepResidue::ExtractResidueCoordinateType(std::istream& ss)
        {
            std::string s;
            ss >> s;
            properties.coordinateType = (s == "XYZ") ? kXYZ : kINT;
        }

        /// Return output format from a stream line which is the 3rd column of the 3rd line in each residue section
        void PrepResidue::ExtractResidueOutputFormat(std::istream& ss)
        {
            int val;
            ss >> val;
            properties.outputFormat = (val == 1) ? kBinary : kFormatted;
        }

        /// Return geometry type from a stream line which is the first column of the 4th line in each residue section
        void PrepResidue::ExtractResidueGeometryType(std::istream& ss)
        {
            std::string s;
            ss >> s;
            properties.geometryType = (s == "CHANGE") ? kGeometryChange : kGeometryCorrect;
        }

        /// Return dummy atom omission from a stream line which is the 2nd column of the 4th line in each residue
        /// section
        void PrepResidue::ExtractResidueDummyAtomOmission(std::istream& ss)
        {
            std::string s;
            ss >> s;
            properties.dummyAtomOmission = (s == "NOMIT") ? kNomit : kOmit;
        }

        /// Return dummy atom position from a stream line which is the 4th column of the 4th line in each residue
        /// section
        void PrepResidue::ExtractResidueDummyAtomPosition(std::istream& ss)
        {
            std::string s;
            ss >> s;
            properties.dummyAtomPosition = (s == "ALL") ? kPositionAll : kPositionBeg;
        }

        /// Return a corresponding title from a stream line which may appear in each residue section
        SectionType PrepResidue::ExtractSectionType(std::string& line)
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
        void PrepResidue::ExtractLoops(std::istream& in_file)
        {
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
                properties.loops.push_back({atom_names[0], atom_names[1]});
            }
            return;
        }

        /// Parse the improper dihedral section of each residue section and return a std::vector of improper dihedrals
        void PrepResidue::ExtractImproperDihedral(std::istream& in_file)
        {
            std::string line;
            std::stringstream ss;
            std::vector<Dihedral> dihedrals;
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
            properties.improperDihedrals = dihedrals;
        }

        //////////////////////////////////////////////////////////
        //                     FUNCTIONS                        //
        //////////////////////////////////////////////////////////
        double PrepResidue::CalculatePrepResidueCharge()
        {
            double residue_charge = 0.0;
            for (auto& atom : this->getAtoms())
            {
                residue_charge += atom->getCharge();
            }
            return residue_charge;
        }

        //////////////////////////////////////////////////////////
        //                     DISPLAY FUNCTIONS                //
        //////////////////////////////////////////////////////////
        std::string PrepResidue::toString() const
        {
            std::stringstream out;
            out << "Title: " << properties.title << std::endl;
            out << std::setw(10) << "ResName" << std::setw(10) << "CrdType" << std::setw(10) << "Output"
                << std::setw(10) << "GeoType" << std::setw(15) << "DummyOmission" << std::setw(12) << "DummyType"
                << std::setw(12) << "DummyPos" << std::setw(10) << "Charge" << std::endl;
            out << std::setw(10) << this->getName();

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
            for (auto& improperDihedral : this->properties.improperDihedrals)
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

        void PrepResidue::Write(std::ostream& stream)
        {
            stream << this->properties.title << std::endl
                   << std::endl
                   << std::left << std::setw(4) << this->getName() << " " << std::right << std::setw(3)
                   << this->GetStringFormatOfCoordinateType(this->properties.coordinateType) << " " << std::setw(1)
                   << this->properties.outputFormat << std::endl
                   << this->GetStringFormatOfGeometryType(this->properties.geometryType) << " "
                   << this->GetStringFormatOfDummyAtomOmission(this->properties.dummyAtomOmission) << " "
                   << this->properties.dummyAtomType << " "
                   << this->GetStringFormatOfDummyAtomPosition(this->properties.dummyAtomPosition) << std::endl
                   << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->properties.charge
                   << std::endl;
            for (auto& atom : this->getAtoms())
            {
                util::erratic_cast<PrepAtom*>(atom)->Write(stream);
            }
            stream << std::endl;
            if (this->properties.improperDihedrals.size() > 0)
            {
                stream << "IMPROPER" << std::endl;
                for (auto& dihedral : this->properties.improperDihedrals)
                {
                    stream << std::left << std::setw(4) << dihedral.at(0) << std::left << std::setw(4) << dihedral.at(1)
                           << std::left << std::setw(4) << dihedral.at(2) << std::left << std::setw(4) << dihedral.at(3)
                           << std::endl;
                }
                stream << std::endl;
            }
            if (this->properties.loops.size() > 0)
            {
                stream << "LOOP" << std::endl;
                for (auto& loop : this->properties.loops)
                {
                    stream << loop.first << " " << loop.second << std::endl;
                }
            }
            stream << std::endl << "DONE" << std::endl;
        }
    } // namespace prep
} // namespace gmml
