#include "include/writers/print.hpp"

#include "include/CentralDataStructure/residue.hpp"
#include "include/readers/Pdb/pdbModel.hpp"

#include <iomanip>
#include <ostream>

namespace gmml
{
    void print(std::ostream& out, const Coordinate& coord)
    {
        out << std::setw(10) << coord.GetX() << ", " << std::setw(10) << coord.GetY() << ", " << std::setw(10)
            << coord.GetZ();
    }

    void print(std::ostream& out, const Atom& atom) { out << atom.getName() << ", "; }

    void print(std::ostream& out, const Residue& residue)
    {
        out << "Name: " << residue.getName() << ", index:" << residue.getIndex() << ":\n";
        for (auto& atom : residue.getAtoms())
        {
            print(out, *atom);
            out << "\n";
        }
        out << "\n";
    }
} // namespace gmml
