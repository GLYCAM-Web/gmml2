#include "include/CentralDataStructure/assembly.hpp"

#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/metadata/atomicBonds.hpp" // bondIfClose
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"

namespace gmml
{
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTORS                   //
    //////////////////////////////////////////////////////////
    Assembly::Assembly(std::vector<Residue*>& residues) : Node<Assembly>(glygraph::invalid, {})
    {
        this->addMolecule(std::make_unique<Molecule>(residues));
    }

    // Move Ctor
    Assembly::Assembly(Assembly&& other) noexcept : glygraph::Node<Assembly>(other)
    {
        molecules_ = std::move(other.molecules_);
        number_ = std::move(other.number_);
    }

    // Copy Ctor
    Assembly::Assembly(const Assembly& other) : glygraph::Node<Assembly>(other), number_(other.number_)
    {
        for (auto& molecule : other.molecules_)
        {
            molecules_.push_back(std::make_unique<Molecule>((*molecule.get())));
        }
    }

    // Move and Copy assignment operator
    Assembly& Assembly::operator=(Assembly other)
    {
        glygraph::Node<Assembly>::operator=(other); // ToDo fuck.
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
            residues.insert(
                residues.end(),
                std::make_move_iterator(currentMoleculeResidues.begin()),
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
            atoms.insert(
                atoms.end(), // Concatenates the vectors. currentResidueAtoms isn't left in a defined state but
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

    const Atom* Assembly::findAtom(uint serialNumber) const
    {
        std::vector<Atom*> atoms = getAtoms();
        std::vector<uint> numbers = atomNumbers(atoms);
        size_t index = util::indexOf(numbers, serialNumber);
        return index < atoms.size() ? atoms[index] : nullptr;
    }
} // namespace gmml
