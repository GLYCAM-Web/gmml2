#include "include/CentralDataStructure/cdsFunctions.hpp"

#include "include/CentralDataStructure/atom.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/elements.hpp"
#include "include/util/containerTypes.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"

#include <array>
#include <vector>

namespace gmml
{
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

    std::vector<Sphere> atomCoordinatesWithRadii(
        const util::SparseVector<double>& elementRadii, const std::vector<Atom*>& atoms)
    {
        std::vector<Sphere> spheres;
        spheres.reserve(atoms.size());
        for (auto& atom : atoms)
        {
            Element element = atom->cachedElement();
            if (!elementRadii.hasValue[element])
            {
                std::string message = "No valid radius for element: " + std::to_string(element);
                util::log(__LINE__, __FILE__, util::ERR, message);
                throw std::runtime_error(message);
            }
            spheres.push_back({elementRadii.values[element], atom->coordinate()});
        }
        return spheres;
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
