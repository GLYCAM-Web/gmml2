#include "includes/CentralDataStructure/Writers/print.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbModel.hpp"

#include <iomanip>
#include <ostream>

void cds::print(std::ostream& out, const cds::Coordinate& coord)
{
    out << std::setw(10) << coord.GetX() << ", " << std::setw(10) << coord.GetY() << ", " << std::setw(10)
        << coord.GetZ();
}

void cds::print(std::ostream& out, const cds::Atom& atom)
{
    out << atom.getName() << ", ";
}

void cds::print(std::ostream& out, const cds::Residue& residue)
{
    out << "Name: " << residue.getName() << ", index:" << residue.getIndex() << ":\n";
    for (auto& atom : residue.getAtoms())
    {
        print(out, *atom);
        out << "\n";
    }
    out << "\n";
}
