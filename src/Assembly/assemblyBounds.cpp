#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"

#include "includes/CodeUtils/containers.hpp"

namespace assembly
{
    void updateResidueBounds(const Graph& graph, Bounds& bounds, size_t index)
    {
        bounds.residues[index] =
            cds::boundingSphere(codeUtils::indicesToValues(bounds.atoms, residueAtoms(graph, index)));
    }

    void updateResidueMoleculeBounds(const Graph& graph, Bounds& bounds, size_t index)
    {
        size_t moleculeIndex = graph.residueMolecule[index];
        bounds.molecules[moleculeIndex] =
            cds::boundingSphereIncluding(bounds.molecules[moleculeIndex], bounds.residues[index]);
    }

    void updateMoleculeBounds(const Graph& graph, Bounds& bounds, size_t index)
    {
        bounds.molecules[index] =
            cds::boundingSphere(codeUtils::indicesToValues(bounds.residues, moleculeResidues(graph, index)));
    }

    void updateBoundsContainingAtoms(const Graph& graph, Bounds& bounds, const std::vector<size_t>& atoms)
    {
        std::vector<bool> updateResidue(graph.residueCount, false);
        std::vector<bool> updateMolecule(graph.moleculeCount, false);
        for (size_t atom : atoms)
        {
            updateResidue[graph.atomResidue[atom]] = true;
        }
        for (size_t n : codeUtils::boolsToIndices(updateResidue))
        {
            updateMolecule[graph.residueMolecule[n]] = true;
            updateResidueBounds(graph, bounds, n);
        }
        for (size_t n : codeUtils::boolsToIndices(updateMolecule))
        {
            updateMoleculeBounds(graph, bounds, n);
        }
    }
} // namespace assembly
