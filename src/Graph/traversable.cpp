#include "includes/Graph/traversable.hpp"
#include "includes/Graph/graphDataLayer.hpp"
#include "includes/Graph/types.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>

namespace graph
{
    void initTraversableStructure(TraversableAssemblyStructure& structure, GraphDataLayer& layer,
                                  graph::Database& database)
    {
        structure.atomData        = {&layer.atomCoordinates, &layer.atomNames, &layer.atomElements, &layer.atomCharges};
        structure.residueData     = {&layer.residueNames};
        structure.moleculeData    = {&layer.moleculeNames};
        structure.atomLinkageData = {&layer.bondTypes};
        structure.residueLinkageData  = {&layer.residueLinkages};
        structure.moleculeLinkageData = {&layer.moleculeLinkages};

        structure.atomGraph     = graph::identity(database);
        structure.residueGraph  = graph::quotient(database, layer.atomResidue);
        structure.moleculeGraph = graph::quotient(graph::asData(structure.residueGraph), layer.residueMolecule);
        Graph& atomGraph        = structure.atomGraph;
        Graph& residueGraph     = structure.residueGraph;
        Graph& moleculeGraph    = structure.moleculeGraph;

        std::vector<size_t>& atomGraphIndices = atomGraph.nodes.indices;
        size_t atomCount                      = atomGraphIndices.size();
        structure.atoms.reserve(atomCount);
        for (size_t n = 0; n < atomCount; n++)
        {
            structure.atoms.push_back({atomGraphIndices[n], &structure.atomData, n, &structure.atomGraphData});
        }

        std::vector<size_t>& residueGraphIndices = residueGraph.nodes.indices;
        size_t residueCount                      = residueGraphIndices.size();
        structure.residues.reserve(residueCount);
        for (size_t n = 0; n < residueCount; n++)
        {
            structure.residues.push_back(
                {residueGraphIndices[n], &structure.residueData, n, &structure.residueGraphData});
        }

        std::vector<size_t>& moleculeGraphIndices = moleculeGraph.nodes.indices;
        size_t moleculeCount                      = moleculeGraphIndices.size();
        structure.molecules.reserve(moleculeCount);
        for (size_t n = 0; n < moleculeCount; n++)
        {
            structure.molecules.push_back(
                {moleculeGraphIndices[n], &structure.moleculeData, n, &structure.moleculeGraphData});
        }

        std::vector<size_t>& atomLinkageGraphIndices = atomGraph.edges.indices;
        size_t atomLinkageCount                      = atomLinkageGraphIndices.size();
        structure.atomLinkages.reserve(atomLinkageCount);
        for (size_t n = 0; n < atomLinkageCount; n++)
        {
            structure.atomLinkages.push_back(
                {atomLinkageGraphIndices[n], &structure.atomLinkageData, n, &structure.atomLinkageGraphData});
        }

        std::vector<size_t>& residueLinkageGraphIndices = residueGraph.edges.indices;
        size_t residueLinkageCount                      = residueLinkageGraphIndices.size();
        structure.residueLinkages.reserve(residueLinkageCount);
        for (size_t n = 0; n < residueLinkageCount; n++)
        {
            structure.residueLinkages.push_back(
                {residueLinkageGraphIndices[n], &structure.residueLinkageData, n, &structure.residueLinkageGraphData});
        }

        std::vector<size_t>& moleculeLinkageGraphIndices = moleculeGraph.edges.indices;
        size_t moleculeLinkageCount                      = moleculeLinkageGraphIndices.size();
        structure.moleculeLinkages.reserve(moleculeLinkageCount);
        for (size_t n = 0; n < moleculeLinkageCount; n++)
        {
            structure.moleculeLinkages.push_back({moleculeLinkageGraphIndices[n], &structure.moleculeLinkageData, n,
                                                  &structure.moleculeLinkageGraphData});
        }

        structure.atomGraphData.residues.reserve(atomCount);
        for (size_t n = 0; n < atomCount; n++)
        {
            structure.atomGraphData.residues.push_back(&structure.residues[layer.atomResidue[n]]);
        }

        structure.residueGraphData.molecules.reserve(residueCount);
        for (size_t n = 0; n < residueCount; n++)
        {
            structure.residueGraphData.molecules.push_back(&structure.molecules[layer.residueMolecule[n]]);
        }

        structure.atomGraphData.neighbors.reserve(atomCount);
        structure.atomGraphData.linkages.reserve(atomCount);
        for (size_t n = 0; n < atomCount; n++)
        {
            std::vector<size_t>& nodes = atomGraph.nodes.nodeAdjacencies[n];
            std::vector<size_t>& edges = atomGraph.nodes.edgeAdjacencies[n];
            size_t adjCount            = edges.size();
            std::vector<Atom*> neighbors;
            neighbors.reserve(adjCount);
            std::vector<AtomLinkage*> linkages;
            linkages.reserve(adjCount);
            for (size_t k = 0; k < adjCount; k++)
            {
                neighbors.push_back(&structure.atoms[nodes[k]]);
                linkages.push_back(&structure.atomLinkages[edges[k]]);
            }
            structure.atomGraphData.neighbors.push_back(neighbors);
            structure.atomGraphData.linkages.push_back(linkages);
        }

        structure.residueGraphData.neighbors.reserve(residueCount);
        structure.residueGraphData.linkages.reserve(residueCount);
        structure.residueGraphData.atoms.reserve(residueCount);
        for (size_t n = 0; n < residueCount; n++)
        {
            std::vector<size_t>& nodes    = residueGraph.nodes.nodeAdjacencies[n];
            std::vector<size_t>& edges    = residueGraph.nodes.edgeAdjacencies[n];
            std::vector<size_t>& elements = residueGraph.nodes.elements[n];
            size_t adjCount               = edges.size();
            std::vector<Residue*> neighbors;
            neighbors.reserve(adjCount);
            std::vector<ResidueLinkage*> linkages;
            linkages.reserve(adjCount);
            std::vector<Atom*> atoms;
            atoms.reserve(elements.size());
            for (size_t k = 0; k < adjCount; k++)
            {
                neighbors.push_back(&structure.residues[nodes[k]]);
                linkages.push_back(&structure.residueLinkages[edges[k]]);
            }
            for (size_t k = 0; k < elements.size(); k++)
            {
                atoms.push_back(&structure.atoms[elements[k]]);
            }
            structure.residueGraphData.neighbors.push_back(neighbors);
            structure.residueGraphData.linkages.push_back(linkages);
            structure.residueGraphData.atoms.push_back(atoms);
        }

        structure.moleculeGraphData.neighbors.reserve(moleculeCount);
        structure.moleculeGraphData.linkages.reserve(moleculeCount);
        structure.moleculeGraphData.atoms.reserve(moleculeCount);
        structure.moleculeGraphData.residues.reserve(moleculeCount);
        for (size_t n = 0; n < moleculeCount; n++)
        {
            std::vector<size_t>& nodes    = moleculeGraph.nodes.nodeAdjacencies[n];
            std::vector<size_t>& edges    = moleculeGraph.nodes.edgeAdjacencies[n];
            std::vector<size_t>& elements = moleculeGraph.nodes.elements[n];
            size_t adjCount               = edges.size();
            std::vector<Molecule*> neighbors;
            neighbors.reserve(adjCount);
            std::vector<MoleculeLinkage*> linkages;
            linkages.reserve(adjCount);
            std::vector<Residue*> residues;
            residues.reserve(elements.size());
            std::vector<Atom*> atoms;
            for (size_t k = 0; k < adjCount; k++)
            {
                neighbors.push_back(&structure.molecules[nodes[k]]);
                linkages.push_back(&structure.moleculeLinkages[edges[k]]);
            }
            for (size_t k = 0; k < elements.size(); k++)
            {
                residues.push_back(&structure.residues[elements[k]]);
                codeUtils::insertInto(atoms, structure.residueGraphData.atoms[elements[k]]);
            }
            structure.moleculeGraphData.neighbors.push_back(neighbors);
            structure.moleculeGraphData.linkages.push_back(linkages);
            structure.moleculeGraphData.atoms.push_back(atoms);
            structure.moleculeGraphData.residues.push_back(residues);
        }

        structure.atomLinkageGraphData.atoms.reserve(atomLinkageCount);
        for (size_t n = 0; n < atomLinkageCount; n++)
        {
            std::array<size_t, 2> nodes = atomGraph.edges.nodeAdjacencies[n];
            structure.atomLinkageGraphData.atoms.push_back({&structure.atoms[nodes[0]], &structure.atoms[nodes[1]]});
        }

        structure.residueLinkageGraphData.residues.reserve(residueLinkageCount);
        for (size_t n = 0; n < residueLinkageCount; n++)
        {
            std::array<size_t, 2> nodes = residueGraph.edges.nodeAdjacencies[n];
            structure.residueLinkageGraphData.residues.push_back(
                {&structure.residues[nodes[0]], &structure.residues[nodes[1]]});
        }

        structure.moleculeLinkageGraphData.molecules.reserve(moleculeLinkageCount);
        for (size_t n = 0; n < moleculeLinkageCount; n++)
        {
            std::array<size_t, 2> nodes = moleculeGraph.edges.nodeAdjacencies[n];
            structure.moleculeLinkageGraphData.molecules.push_back(
                {&structure.molecules[nodes[0]], &structure.molecules[nodes[1]]});
        }
    }
} // namespace graph
