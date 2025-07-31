#include "include/CentralDataStructure/cdsFunctions.hpp"

#include "include/CentralDataStructure/assembly.hpp"
#include "include/CentralDataStructure/atom.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/elements.hpp"
#include "include/pdb/pdbData.hpp"
#include "include/pdb/pdbResidue.hpp"
#include "include/pdb/pdbResidueId.hpp"
#include "include/util/containerTypes.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"

#include <array>
#include <vector>

namespace gmml
{
    std::vector<Atom*> getAtoms(const std::vector<Assembly*>& assemblies)
    {
        std::vector<Atom*> result;
        for (auto& assembly : assemblies)
        {
            util::insertInto(result, assembly->getAtoms());
        }
        return result;
    }

    std::vector<Residue*> getResidues(const std::vector<Assembly*>& assemblies)
    {
        std::vector<Residue*> result;
        for (auto& assembly : assemblies)
        {
            util::insertInto(result, assembly->getResidues());
        }
        return result;
    }

    size_t residueSelector(const pdb::PdbData& data, const pdb::ResidueId& residueId, uint modelNumber)
    {
        size_t n = util::indexOf(data.assemblies.numbers, modelNumber);
        if (n < data.assemblies.numbers.size())
        {
            return residueSelector(data, assemblyResidues(data.indices, n), residueId);
        }
        else
        {
            return data.indices.residueCount;
        }
    }

    size_t residueSelector(const pdb::PdbData& data, std::vector<size_t> residueIds, const pdb::ResidueId& queryId)
    { // I'm using empty() to mean that it could be anything.
        auto eq = [](const std::string& str, const std::string& comp) { return str.empty() || str == comp; };

        for (size_t residueId : residueIds)
        {
            pdb::ResidueId id = pdb::pdbResidueId(data, residueId);
            if (eq(queryId.residueName, id.residueName) && eq(queryId.sequenceNumber, id.sequenceNumber) &&
                eq(queryId.insertionCode, id.insertionCode) && eq(queryId.chainId, id.chainId))
            {
                return residueId;
            }
        }
        return data.indices.residueCount;
    }

    std::vector<uint> serializedNumberVector(const std::vector<bool>& included)
    {
        std::vector<uint> result;
        result.reserve(included.size());
        uint count = 0;
        for (bool b : included)
        {
            result.push_back(b ? count + 1 : 0);
            count += b;
        }
        return result;
    }

    std::vector<uint> serializedNumberVector(size_t count)
    {
        return serializedNumberVector(std::vector<bool>(count, true));
    }

    size_t atomVectorIndex(const std::vector<Atom*>& atoms, Atom* find)
    {
        size_t index = util::indexOf(atoms, find);
        if (index == atoms.size())
        {
            throw std::runtime_error("atom missing from data in off writer data");
        }
        // index equals offset from start of vector
        return index;
    }

    std::vector<uint> atomNumbers(const std::vector<Atom*>& atoms)
    {
        std::vector<uint> result;
        result.reserve(atoms.size());
        for (auto& atom : atoms)
        {
            result.push_back(atom->getNumber());
        }
        return result;
    }

    std::vector<std::string> atomNames(const std::vector<Atom*>& atoms)
    {
        std::vector<std::string> result;
        result.reserve(atoms.size());
        for (auto& atom : atoms)
        {
            result.push_back(atom->getName());
        }
        return result;
    }

    std::vector<std::string> atomElementStrings(const std::vector<Atom*>& atoms)
    {
        std::vector<std::string> result;
        result.reserve(atoms.size());
        for (auto& atom : atoms)
        {
            result.push_back(atom->getElement());
        }
        return result;
    }

    std::vector<Element> atomElements(const std::vector<Atom*>& atoms)
    {
        std::vector<Element> result;
        result.reserve(atoms.size());
        for (auto& atom : atoms)
        {
            result.push_back(atom->cachedElement());
        }
        return result;
    }

    std::vector<uint> atomAtomicNumbers(const std::vector<Atom*>& atoms)
    {
        std::vector<uint> result;
        result.reserve(atoms.size());
        for (auto& atom : atoms)
        {
            result.push_back(atom->getAtomicNumber());
        }
        return result;
    }

    std::vector<Coordinate> atomCoordinates(const std::vector<Atom*>& atoms)
    {
        std::vector<Coordinate> coordinates;
        coordinates.reserve(atoms.size());
        for (auto& atom : atoms)
        {
            coordinates.push_back(atom->coordinate());
        }
        return coordinates;
    }

    void setAtomCoordinates(std::vector<Atom*>& atoms, const std::vector<Coordinate>& coordinates)
    {
        if (atoms.size() != coordinates.size())
        {
            throw std::runtime_error("Trying to set atom vector coordinates with wrong sized coordinate vector");
        }
        for (size_t n = 0; n < atoms.size(); n++)
        {
            atoms[n]->setCoordinate(coordinates[n]);
        }
    }

    std::vector<std::string> atomTypes(const std::vector<Atom*>& atoms)
    {
        std::vector<std::string> result;
        result.reserve(atoms.size());
        for (auto& atom : atoms)
        {
            result.push_back(atom->getType());
        }
        return result;
    }

    std::vector<double> atomCharges(const std::vector<Atom*>& atoms)
    {
        std::vector<double> result;
        result.reserve(atoms.size());
        for (auto& atom : atoms)
        {
            result.push_back(atom->getCharge());
        }
        return result;
    }

    std::vector<bool> atomVisibility(const std::vector<Atom*>& atoms)
    {
        std::vector<bool> result;
        result.reserve(atoms.size());
        for (auto& atom : atoms)
        {
            result.push_back(atom->isVisible());
        }
        return result;
    }

    std::vector<uint> residueNumbers(const std::vector<Residue*>& residues)
    {
        std::vector<uint> result;
        result.reserve(residues.size());
        for (auto& residue : residues)
        {
            result.push_back(residue->getNumber());
        }
        return result;
    }

    std::vector<std::string> residueNames(const std::vector<Residue*>& residues)
    {
        std::vector<std::string> result;
        result.reserve(residues.size());
        for (auto& residue : residues)
        {
            result.push_back(residue->getName());
        }
        return result;
    }

    std::vector<std::string> residueStringIds(const std::vector<Residue*>& residues)
    {
        std::vector<std::string> result;
        result.reserve(residues.size());
        for (auto& residue : residues)
        {
            result.push_back(residueStringId(residue));
        }
        return result;
    }

    std::string truncatedResidueName(const Residue* residue) { return residue->getName().substr(0, 3); }

    std::vector<std::string> truncatedResidueNames(const std::vector<Residue*>& residues)
    {
        std::vector<std::string> result;
        result.reserve(residues.size());
        for (auto& residue : residues)
        {
            result.push_back(truncatedResidueName(residue));
        }
        return result;
    }

    std::vector<ResidueType> residueTypes(const std::vector<Residue*>& residues)
    {
        std::vector<ResidueType> result;
        result.reserve(residues.size());
        for (auto& residue : residues)
        {
            result.push_back(residue->GetType());
        }
        return result;
    }

    void addBond(Atom* atom, Atom* otherAtom) { atom->addNeighbor("atomicBond", otherAtom); }
} // namespace gmml
