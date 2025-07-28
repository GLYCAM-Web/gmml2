#ifndef INCLUDE_CENTRALDATASTRUCTURE_CDSFUNCTIONS_HPP
#define INCLUDE_CENTRALDATASTRUCTURE_CDSFUNCTIONS_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/elements.hpp"
#include "include/metadata/residueTypes.hpp"
#include "include/pdb/pdbData.hpp"
#include "include/pdb/pdbResidueId.hpp"
#include "include/util/containerTypes.hpp"

#include <array>
#include <vector>

namespace gmml
{
    template<typename T> void serializeNumbers(std::vector<T*> elements)
    {
        for (uint n = 0; n < elements.size(); n++)
        {
            elements[n]->setNumber(n + 1);
        }
        return;
    }

    std::vector<Atom*> getAtoms(const std::vector<Assembly*>& assemblies);
    std::vector<Residue*> getResidues(const std::vector<Assembly*>& assemblies);
    size_t residueSelector(const pdb::PdbData& data, const pdb::ResidueId& residueId, const int modelNumber = 0);
    size_t residueSelector(const pdb::PdbData& data, std::vector<size_t> residueIds, const pdb::ResidueId& queryId);
    std::vector<uint> serializedNumberVector(const std::vector<bool>& included);
    std::vector<uint> serializedNumberVector(size_t count);
    size_t atomVectorIndex(const std::vector<Atom*>& atoms, Atom* find);
    std::vector<uint> atomNumbers(const std::vector<Atom*>& atoms);
    std::vector<std::string> atomNames(const std::vector<Atom*>& atoms);
    std::vector<std::string> atomElementStrings(const std::vector<Atom*>& atoms);
    std::vector<Element> atomElements(const std::vector<Atom*>& atoms);
    std::vector<uint> atomAtomicNumbers(const std::vector<Atom*>& atoms);
    std::vector<Coordinate> atomCoordinates(const std::vector<Atom*>& atoms);
    void setAtomCoordinates(std::vector<Atom*>& atoms, const std::vector<Coordinate>& coordinates);
    std::vector<Sphere> atomCoordinatesWithRadii(
        const util::SparseVector<double>& elementRadii, const std::vector<Atom*>& atoms);
    std::vector<std::string> atomTypes(const std::vector<Atom*>& atoms);
    std::vector<double> atomCharges(const std::vector<Atom*>& atoms);
    std::vector<bool> atomVisibility(const std::vector<Atom*>& atoms);
    std::vector<uint> residueNumbers(const std::vector<Residue*>& residues);
    std::vector<std::string> residueNames(const std::vector<Residue*>& residues);
    std::vector<std::string> residueStringIds(const std::vector<Residue*>& residues);
    std::string truncatedResidueName(const Residue* residue);
    std::vector<std::string> truncatedResidueNames(const std::vector<Residue*>& residues);
    std::vector<ResidueType> residueTypes(const std::vector<Residue*>& residues);
    void addBond(Atom* atom, Atom* otherAtom);
} // namespace gmml

#endif
