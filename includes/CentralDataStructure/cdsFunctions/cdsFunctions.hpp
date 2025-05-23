#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CodeUtils/containerTypes.hpp"
#include "includes/CodeUtils/references.hpp"

#include <array>
#include <vector>

namespace cds
{
    template<typename T> void serializeNumbers(std::vector<T*> elements)
    {
        for (size_t n = 0; n < elements.size(); n++)
        {
            elements[n]->setNumber(n + 1);
        }
        return;
    }

    std::vector<int> serializedNumberVector(const std::vector<bool>& included);
    std::vector<int> serializedNumberVector(size_t count);
    size_t atomVectorIndex(const std::vector<cds::Atom*>& atoms, cds::Atom* find);
    std::vector<int> atomNumbers(const std::vector<Atom*>& atoms);
    std::vector<std::string> atomNames(const std::vector<Atom*>& atoms);
    std::vector<std::string> atomElementStrings(const std::vector<Atom*>& atoms);
    std::vector<MolecularMetadata::Element> atomElements(const std::vector<Atom*>& atoms);
    std::vector<int> atomAtomicNumbers(const std::vector<Atom*>& atoms);
    std::vector<Coordinate> atomCoordinates(const std::vector<Atom*>& atoms);
    std::vector<Sphere> atomCoordinatesWithRadii(const codeUtils::SparseVector<double>& elementRadii,
                                                 const std::vector<Atom*>& atoms);
    std::vector<CoordinateReference> atomCoordinateReferences(std::vector<Atom*>& atoms);
    std::vector<std::string> atomTypes(const std::vector<Atom*>& atoms);
    std::vector<double> atomCharges(const std::vector<Atom*>& atoms);
    std::vector<bool> atomVisibility(const std::vector<Atom*>& atoms);
    std::vector<int> residueNumbers(const std::vector<Residue*>& residues);
    std::vector<std::string> residueNames(const std::vector<Residue*>& residues);
    std::vector<std::string> residueStringIds(const std::vector<Residue*>& residues);
    std::string truncatedResidueName(const Residue* residue);
    std::vector<std::string> truncatedResidueNames(const std::vector<Residue*>& residues);
    std::vector<ResidueType> residueTypes(const std::vector<Residue*>& residues);
} // namespace cds
#endif
