#ifndef INCLUDES_CENTRALDATASTRUCTURE_GRAPH_TRAVERSABLE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GRAPH_TRAVERSABLE_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/Graph/types.hpp"
#include "includes/Graph/manipulation.hpp"
#include "includes/Graph/graphDataLayer.hpp"

#include <cstddef>
#include <array>
#include <optional>
#include <vector>
#include <string>

namespace graph
{
    class Atom;
    class AtomLinkage;
    class Residue;
    class ResidueLinkage;
    class Molecule;
    class MoleculeLinkage;

    struct AtomData
    {
        std::vector<cds::Coordinate>* coordinates;
        std::vector<std::string>* names;
        std::vector<std::string>* elements;
        std::vector<double>* charges;
    };

    struct AtomGraphData
    {
        std::vector<std::vector<AtomLinkage*>> linkages;
        std::vector<std::vector<Atom*>> neighbors;
        std::vector<Residue*> residues;
    };

    struct AtomLinkageData
    {
        std::vector<BondType>* bondTypes;
    };

    struct AtomLinkageGraphData
    {
        std::vector<std::array<Atom*, 2>> atoms;
    };

    struct ResidueData
    {
        std::vector<std::string>* names;
    };

    struct ResidueGraphData
    {
        std::vector<std::vector<ResidueLinkage*>> linkages;
        std::vector<std::vector<Residue*>> neighbors;
        std::vector<std::vector<Atom*>> atoms;
        std::vector<Molecule*> molecules;
    };

    struct ResidueLinkageData
    {
        std::vector<std::optional<ResidueLinkageStruct>>* structs;
    };

    struct ResidueLinkageGraphData
    {
        std::vector<std::array<Residue*, 2>> residues;
    };

    struct MoleculeData
    {
        std::vector<std::string>* names;
    };

    struct MoleculeGraphData
    {
        std::vector<std::vector<MoleculeLinkage*>> linkages;
        std::vector<std::vector<Molecule*>> neighbors;
        std::vector<std::vector<Atom*>> atoms;
        std::vector<std::vector<Residue*>> residues;
    };

    struct MoleculeLinkageData
    {
        std::vector<std::optional<MoleculeLinkageStruct>>* structs;
    };

    struct MoleculeLinkageGraphData
    {
        std::vector<std::array<Molecule*, 2>> molecules;
    };

    class AtomLinkage
    {
      public:
        AtomLinkage(size_t dataIndex_, AtomLinkageData* data_, size_t graphIndex_, AtomLinkageGraphData* graph_)
            : dataIndex(dataIndex_), data(data_), graphIndex(graphIndex_), graph(graph_) {};

        inline BondType bondType() const
        {
            return data->bondTypes->at(dataIndex);
        }

        inline std::array<Atom*, 2> atoms() const
        {
            return graph->atoms.at(graphIndex);
        }

      private:
        size_t dataIndex;
        AtomLinkageData* data;
        size_t graphIndex;
        AtomLinkageGraphData* graph;
    };

    class ResidueLinkage
    {
      public:
        ResidueLinkage(size_t dataIndex_, ResidueLinkageData* data_, size_t graphIndex_,
                       ResidueLinkageGraphData* graph_)
            : dataIndex(dataIndex_), data(data_), graphIndex(graphIndex_), graph(graph_) {};

        const ResidueLinkageStruct& stuff() const
        {
            return data->structs->at(dataIndex).value();
        }

        inline Residue* source() const
        {
            return graph->residues.at(graphIndex)[0];
        }

        inline Residue* target() const
        {
            return graph->residues.at(graphIndex)[1];
        }

      private:
        size_t dataIndex;
        ResidueLinkageData* data;
        size_t graphIndex;
        ResidueLinkageGraphData* graph;
    };

    class MoleculeLinkage
    {
      public:
        MoleculeLinkage(size_t dataIndex_, MoleculeLinkageData* data_, size_t graphIndex_,
                        MoleculeLinkageGraphData* graph_)
            : dataIndex(dataIndex_), data(data_), graphIndex(graphIndex_), graph(graph_) {};

        const MoleculeLinkageStruct& stuff() const
        {
            return data->structs->at(dataIndex).value();
        }

        inline Molecule* source() const
        {
            return graph->molecules.at(graphIndex)[0];
        }

        inline Molecule* target() const
        {
            return graph->molecules.at(graphIndex)[1];
        }

      private:
        size_t dataIndex;
        MoleculeLinkageData* data;
        size_t graphIndex;
        MoleculeLinkageGraphData* graph;
    };

    class Atom
    {
      public:
        Atom(size_t dataIndex_, AtomData* data_, size_t graphIndex_, AtomGraphData* graph_)
            : dataIndex(dataIndex_), data(data_), graphIndex(graphIndex_), graph(graph_) {};

        inline size_t index() const
        {
            return dataIndex;
        }

        inline cds::Coordinate coordinate() const
        {
            return data->coordinates->at(dataIndex);
        }

        inline void setCoordinate(const cds::Coordinate& coord)
        {
            data->coordinates->at(dataIndex) = coord;
        }

        inline std::string& name() const
        {
            return data->names->at(dataIndex);
        }

        inline const std::string& element() const
        {
            return data->elements->at(dataIndex);
        }

        inline double charge() const
        {
            return data->charges->at(dataIndex);
        }

        inline const std::vector<AtomLinkage*>& linkages() const
        {
            return graph->linkages.at(graphIndex);
        }

        inline const std::vector<Atom*>& neighbors() const
        {
            return graph->neighbors.at(graphIndex);
        }

        inline const Residue* residue() const
        {
            return graph->residues.at(graphIndex);
        }

      private:
        size_t dataIndex;
        AtomData* data;
        size_t graphIndex;
        AtomGraphData* graph;
    };

    class Residue
    {
      public:
        Residue(size_t dataIndex_, ResidueData* data_, size_t graphIndex_, ResidueGraphData* graph_)
            : dataIndex(dataIndex_), data(data_), graphIndex(graphIndex_), graph(graph_) {};

        inline size_t index() const
        {
            return dataIndex;
        }

        inline std::string& name() const
        {
            return data->names->at(dataIndex);
        }

        inline const std::vector<ResidueLinkage*>& linkages() const
        {
            return graph->linkages.at(graphIndex);
        }

        inline const std::vector<Residue*>& neighbors() const
        {
            return graph->neighbors.at(graphIndex);
        }

        inline const std::vector<Atom*>& atoms() const
        {
            return graph->atoms.at(graphIndex);
        }

        inline Molecule* molecule() const
        {
            return graph->molecules.at(graphIndex);
        }

      private:
        size_t dataIndex;
        ResidueData* data;
        size_t graphIndex;
        ResidueGraphData* graph;
    };

    class Molecule
    {
      public:
        Molecule(size_t dataIndex_, MoleculeData* data_, size_t graphIndex_, MoleculeGraphData* graph_)
            : dataIndex(dataIndex_), data(data_), graphIndex(graphIndex_), graph(graph_) {};

        inline size_t index() const
        {
            return dataIndex;
        }

        inline std::string name() const
        {
            return data->names->at(dataIndex);
        }

        inline const std::vector<MoleculeLinkage*>& linkages() const
        {
            return graph->linkages.at(graphIndex);
        }

        inline const std::vector<Molecule*>& neighbors() const
        {
            return graph->neighbors.at(graphIndex);
        }

        inline const std::vector<Atom*>& atoms() const
        {
            return graph->atoms.at(graphIndex);
        }

        inline const std::vector<Residue*>& residues() const
        {
            return graph->residues.at(graphIndex);
        }

      private:
        size_t dataIndex;
        MoleculeData* data;
        size_t graphIndex;
        MoleculeGraphData* graph;
    };

    struct TraversableAssemblyStructure
    {
        AtomData atomData;
        ResidueData residueData;
        MoleculeData moleculeData;
        AtomLinkageData atomLinkageData;
        ResidueLinkageData residueLinkageData;
        MoleculeLinkageData moleculeLinkageData;

        graph::Graph atomGraph;
        graph::Graph residueGraph;
        graph::Graph moleculeGraph;

        std::vector<Atom> atoms;
        std::vector<Residue> residues;
        std::vector<Molecule> molecules;
        std::vector<AtomLinkage> atomLinkages;
        std::vector<ResidueLinkage> residueLinkages;
        std::vector<MoleculeLinkage> moleculeLinkages;

        AtomGraphData atomGraphData;
        AtomLinkageGraphData atomLinkageGraphData;
        ResidueGraphData residueGraphData;
        ResidueLinkageGraphData residueLinkageGraphData;
        MoleculeGraphData moleculeGraphData;
        MoleculeLinkageGraphData moleculeLinkageGraphData;
    };

    void initTraversableStructure(TraversableAssemblyStructure& structure, GraphDataLayer& layer, Database& database);
} // namespace graph

#endif
