#ifndef INCLUDES_CENTRALDATASTRUCTURE_MOLECULE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_MOLECULE_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"

#include <vector>
#include <memory> // unique_ptr

namespace cds
{
    class Molecule : public glygraph::Node<Molecule>
    {
      public:
        //////////////////////////////////////////////////////////
        //                    CONSTRUCTOR                       //
        //////////////////////////////////////////////////////////
        Molecule() : Node<Molecule>("cdsMoleculeDefault", {}), number_(0) {};

        Molecule(std::vector<Residue*>& residues);
        Molecule(Molecule&& other) noexcept; // Move Ctor
        Molecule(const Molecule& other);     // Copy Ctor
        Molecule& operator=(Molecule other); // Move and Copy assignment operator

        virtual ~Molecule()
        {} // std::cout << "Molecule default dtor for " << this->getName() << "\n";}

        friend void swap(Molecule& lhs,
                         Molecule& rhs) // ToDo figure out how to put this in cpp file once everything is working.
        {
            using std::swap;
            swap(lhs.residues_, rhs.residues_);
            swap(lhs.number_, rhs.number_); // @suppress("Invalid arguments")
        }

        //////////////////////////////////////////////////////////
        //                    ACCESSOR                          //
        //////////////////////////////////////////////////////////
        inline const int& getNumber()
        {
            return number_;
        }

        std::vector<Atom*> getAtoms() const;
        std::vector<Atom*> mutableAtoms();
        std::vector<Residue*> getResidues() const;

        //////////////////////////////////////////////////////////
        //                    MUTATOR                           //
        //////////////////////////////////////////////////////////
        inline void setNumber(const int i)
        {
            number_ = i;
        }

        void swapResiduePosition(Residue* queryResidue, size_t newPosition);
        //////////////////////////////////////////////////////////
        //                    FUNCTIONS                         //
        //////////////////////////////////////////////////////////
        Residue* addResidue(std::unique_ptr<Residue> myResidue);
        Residue* insertNewResidue(std::unique_ptr<Residue> myResidue, size_t position);
        Residue* insertNewResidue(std::unique_ptr<Residue> myResidue, const Residue& positionReferenceResidue);
        std::vector<std::unique_ptr<Residue>>::iterator findPositionOfResidue(const Residue* queryResidue);
        std::vector<Residue*> getResidues(const std::vector<std::string>& queryNames);
        Residue* getResidue(const std::string& queryName);
        void deleteResidue(Residue*);

      private:
        //////////////////////////////////////////////////////////
        //                    FUNCTIONS                         //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                    ATTRIBUTES                        //
        //////////////////////////////////////////////////////////
        std::vector<std::unique_ptr<Residue>> residues_;
        int number_;
    };
} // namespace cds
#endif
