#ifndef INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP

#include "includes/CentralDataStructure/residueTypes.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp" // getId
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include "includes/CodeUtils/constants.hpp" // iNotSet

#include <vector>
#include <memory>    // unique_ptr
#include <algorithm> //swap
#include <functional>

using cds::Coordinate;

namespace cds
{
    class Residue : public glygraph::Node<Residue>
    {
      public:
        //////////////////////////////////////////////////////////
        //                    CONSTRUCTOR                       //
        //////////////////////////////////////////////////////////
        Residue() : Node<Residue>(glygraph::invalid, {}) {};

        Residue(const std::string& residueName, const Residue* referenceResidue);
        Residue(Residue&& other) noexcept; // Move Ctor
        Residue(const Residue& other);     // Copy Ctor
        Residue& operator=(Residue other); // Move & Copy assignment operator.

        virtual ~Residue()
        {} // std::cout << "Residue dtor for " << this->getName() << std::endl;}// Dtor, virtual so that derived classes
           // dtors will get triggered if possible.

        //////////////////////////////////////////////////////////
        //                    ACCESSOR                          //
        //////////////////////////////////////////////////////////
        inline virtual const std::string& getName() const
        {
            return name_;
        }

        const std::string GetParmName() const;
        std::vector<Atom*> getAtoms() const;
        std::vector<Atom*> getHydrogenAtoms() const;
        std::vector<Atom*> getNonHydrogenAtoms() const;
        std::vector<Atom*> mutableAtoms();
        std::vector<std::string> getAtomNames() const;
        std::string getStringId(std::string moleculeNumber = constants::sNotSet) const;

        inline ResidueType GetType() const
        {
            return type_;
        }

        inline int getNumber() const
        {
            return number_;
        }

        virtual pdb::ResidueId getId() const;

        //////////////////////////////////////////////////////////
        //                    MUTATOR                           //
        //////////////////////////////////////////////////////////
        inline void setName(const std::string& s)
        {
            name_ = s;
        }

        inline void setAtoms(std::vector<std::unique_ptr<Atom>> v)
        {
            atoms_ = std::move(v);
        }

        Atom* addAtom(std::unique_ptr<Atom> myAtom);
        Atom* addAtomToFront(std::unique_ptr<Atom> myAtom);
        bool moveAtomToLastPosition(const Atom* atom);
        bool deleteAtom(const Atom* atom);

        inline void SetType(ResidueType type)
        {
            type_ = type;
        }

        inline void setNumber(int i)
        {
            number_ = i;
        }

        inline size_t atomCount() const
        {
            return atoms_.size();
        }

        //////////////////////////////////////////////////////////
        //                    FUNCTIONS                         //
        //////////////////////////////////////////////////////////
        typename std::vector<std::unique_ptr<Atom>>::iterator FindPositionOfAtom(const Atom* queryAtom);
        Atom* FindAtom(const std::string queryName) const;
        Atom* FindAtom(const int& queryNumber) const;
        bool contains(const Atom* queryAtom) const;
        void MakeDeoxy(const std::string oxygenNumber);
        ResidueType determineType(const std::string& residueName);

        //////////////////////////////////////////////////////////
        //                  OPERATOR OVERLOADING                //
        //////////////////////////////////////////////////////////
        virtual bool operator==(const Residue& rhs) const
        {
            return (this->getIndex() == rhs.getIndex());
        }

        virtual bool operator!=(const Residue& rhs) const
        {
            return (this->getIndex() != rhs.getIndex());
        }

        virtual bool operator>(const Residue& rhs) const
        {
            return (this->getNumber() > rhs.getNumber());
        }

        virtual bool operator<(const Residue& rhs) const
        {
            return (this->getNumber() < rhs.getNumber());
        }

        //////////////////////////////////////////////////////////
        //               Copy-Swap                    //
        //////////////////////////////////////////////////////////
        friend void swap(Residue& lhs, Residue& rhs)
        {
            using std::swap;
            swap(lhs.atoms_, rhs.atoms_);
            swap(lhs.name_, rhs.name_);
            swap(lhs.type_, rhs.type_);
            swap(lhs.number_, rhs.number_);
        }

      private:
        std::vector<Atom*> getAtomsConditional(std::function<bool(Atom*)> condition) const;
        //////////////////////////////////////////////////////////
        //                    ATTRIBUTES                        //
        //////////////////////////////////////////////////////////
        std::vector<std::unique_ptr<Atom>> atoms_;
        std::string name_ = "   ";
        ResidueType type_ = Undefined; // enum Type. See enum above.
        int number_ =
            1; // constants::iNotSet; ToDo: For prep residues a default 1 value is good. Is there a reason not to?
    };
} // namespace cds
#endif
