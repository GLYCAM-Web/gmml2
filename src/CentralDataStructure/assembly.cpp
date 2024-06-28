#include "includes/CentralDataStructure/assembly.hpp"

#include "../../includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "../../includes/CentralDataStructure/Selections/templatedSelections.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/numbers.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp" // bondIfClose

using cds::Assembly;
using cds::Atom;
using cds::Molecule;
using cds::Residue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTORS                   //
//////////////////////////////////////////////////////////
Assembly::Assembly(std::vector<Residue*>& residues) : Node<Assembly>(glygraph::invalid, {})
{
    this->addMolecule(std::make_unique<Molecule>(residues));
}

// Move Ctor
Assembly::Assembly(Assembly&& other) noexcept : glygraph::Node<cds::Assembly>(other)
{
    molecules_ = std::move(other.molecules_);
    number_    = std::move(other.number_);
}

// Copy Ctor
Assembly::Assembly(const Assembly& other) : glygraph::Node<cds::Assembly>(other), number_(other.number_)
{
    for (auto& molecule : other.molecules_)
    {
        molecules_.push_back(std::make_unique<Molecule>((*molecule.get())));
    }
}

// Move and Copy assignment operator
Assembly& Assembly::operator=(Assembly other)
{
    glygraph::Node<cds::Assembly>::operator=(other); // ToDo fuck.
    swap(*this, other);
    return *this;
}

//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
std::vector<Molecule*> Assembly::getMolecules() const
{
    std::vector<Molecule*> molecules;
    molecules.reserve(molecules_.size());
    for (auto& molPtr : molecules_)
    {
        molecules.push_back(molPtr.get());
    }
    return molecules;
}

std::vector<Residue*> Assembly::getResidues() const
{
    std::vector<Residue*> residues;
    for (auto& molPtr : this->getMolecules())
    {
        std::vector<Residue*> currentMoleculeResidues = molPtr->getResidues();
        residues.insert(residues.end(), std::make_move_iterator(currentMoleculeResidues.begin()),
                        std::make_move_iterator(currentMoleculeResidues.end()));
    }
    return residues;
}

std::vector<Atom*> Assembly::getAtoms() const
{
    std::vector<Atom*> atoms;
    for (auto& residue : this->getResidues())
    {
        std::vector<Atom*> currentResidueAtoms = residue->getAtoms();
        atoms.insert(atoms.end(), // Concatenates the vectors. currentResidueAtoms isn't left in a defined state but
                                  // that's ok here.
                     std::make_move_iterator(currentResidueAtoms.begin()),
                     std::make_move_iterator(currentResidueAtoms.end()));
    }
    return atoms;
}

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
Molecule* Assembly::addMolecule(std::unique_ptr<Molecule> myMolecule)
{
    molecules_.push_back(std::move(myMolecule));
    return molecules_.back().get();
}

const Atom* Assembly::findAtom(int serialNumber) const
{
    return codeUtils::findElementWithNumber(this->getAtoms(), serialNumber);
}
