#ifndef INCLUDES_CENTRALDATASTRUCTURE_GRAPH_CDSTOGRAPH_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GRAPH_CDSTOGRAPH_HPP

#include "includes/Graph/graphDataLayer.hpp"
#include "includes/CentralDataStructure/assembly.hpp"

namespace graph
{
    GraphDataLayer toGraphDataLayer(cds::Assembly& assembly);
} // namespace graph

#endif
