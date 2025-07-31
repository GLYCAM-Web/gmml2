#include "include/assembly/assemblyBounds.hpp"

#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/elements.hpp"
#include "include/util/containerTypes.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace gmml
{
    namespace assembly
    {
        void updateResidueBounds(const Graph& graph, Bounds& bounds, size_t index)
        {
            bounds.residues[index] = boundingSphere(util::indicesToValues(bounds.atoms, residueAtoms(graph, index)));
        }

        void updateResidueMoleculeBounds(const Graph& graph, Bounds& bounds, size_t index)
        {
            size_t moleculeIndex = graph.indices.residueMolecule[index];
            bounds.molecules[moleculeIndex] =
                boundingSphereIncluding(bounds.molecules[moleculeIndex], bounds.residues[index]);
        }

        void updateMoleculeBounds(const Graph& graph, Bounds& bounds, size_t index)
        {
            bounds.molecules[index] =
                boundingSphere(util::indicesToValues(bounds.residues, moleculeResidues(graph, index)));
        }

        void updateBoundsContainingAtoms(const Graph& graph, Bounds& bounds, const std::vector<size_t>& selectedAtoms)
        {
            const Indices& indices = graph.indices;
            std::vector<bool> updateResidue(indices.residueCount, false);
            std::vector<bool> updateMolecule(indices.moleculeCount, false);
            for (size_t atom : selectedAtoms)
            {
                updateResidue[indices.atomResidue[atom]] = true;
            }
            for (size_t n = 0; n < indices.residueCount; n++)
            {
                if (updateResidue[n])
                {
                    updateMolecule[indices.residueMolecule[n]] = true;
                    updateResidueBounds(graph, bounds, n);
                }
            }
            for (size_t n = 0; n < indices.moleculeCount; n++)
            {
                if (updateMolecule[n])
                {
                    updateMoleculeBounds(graph, bounds, n);
                }
            }
        }

        std::vector<Sphere> toAtomBounds(
            const util::SparseVector<double>& elementRadii,
            const std::vector<Element>& elements,
            const std::vector<Coordinate>& coordinates)
        {
            std::vector<Sphere> spheres;
            spheres.reserve(elements.size());
            for (size_t n = 0; n < elements.size(); n++)
            {
                Element element = elements[n];
                if (!elementRadii.hasValue[element])
                {
                    std::string message = "No valid radius for element: " + elementName(element);
                    util::log(__LINE__, __FILE__, util::ERR, message);
                    throw std::runtime_error(message);
                }
                spheres.push_back({elementRadii.values[element], coordinates[n]});
            }
            return spheres;
        }

        Bounds toAssemblyBounds(const Graph& graph, const std::vector<Sphere>& atomBounds)
        {
            size_t residueCount = graph.indices.residueCount;
            std::vector<Sphere> residueBounds;
            residueBounds.reserve(residueCount);
            for (size_t n = 0; n < residueCount; n++)
            {
                residueBounds.push_back(boundingSphere(util::indicesToValues(atomBounds, residueAtoms(graph, n))));
            }
            size_t moleculeCount = graph.indices.moleculeCount;
            std::vector<Sphere> moleculeBounds;
            moleculeBounds.reserve(moleculeCount);
            for (size_t n = 0; n < moleculeCount; n++)
            {
                moleculeBounds.push_back(
                    boundingSphere(util::indicesToValues(residueBounds, moleculeResidues(graph, n))));
            }
            return {atomBounds, residueBounds, moleculeBounds};
        }
    } // namespace assembly
} // namespace gmml
