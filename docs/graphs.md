## Graphs
### Data
The graph is key to representing molecular structures. The basic data of a graph is:
* a list of nodes, which are simply their own index
* a list of edges, represented by a pair of indices of the nodes they join

This basic form has the advantage of being easy to iterate over and easy to modify. Adding nodes or edges is as simple as adding an entry to the list. Removing nodes or edges is done by flagging them as removed. Since the graph is index-based, actually removing them from the list would invalidate any index referring to the graph, such as the edge data itself.

### Traversal
The downside of the basic form of graph is that it’s poorly suited for traversal. In order to find the nodes adjacent to a given node, we need to look through the entire list of edges. To remedy this, we can generate traversible representations of the basic graph data. The traversible form extends the basic form with a couple of new columns:

* a list of adjacent nodes, per node
* a list of edges, per node, through which the above are connected

This lets us access and iterate over nodes and edges adjacent to any node, given its index.

The downside of the traversible form, in turn, is that updating it risks being error-prone since information is duplicated in multiple places. For this reason, we’ve chosen to keep the traversible graphs locked from modification. If we need to modify a traversible graph, we instead modify the basic data and recreate the traversible form anew.

### Associated data
The graph contains no data apart from its connectivity. In order to associate data with nodes and edges in the graph, we make use of parallel data structures (vectors, structs of vectors) with a matching number of entries to that of the nodes or edges in the graph, whichever they are associated with. This way, the same index used for iterating over or traversing the graph can be used to access associated data.

### Subgraphs
A subgraph is a traversible graph containing a subset of the nodes and edges of another graph. The node indices in our subgraphs are preserved from the original graph, in order to make data access convenient. The edge indices are not preserved, and for this reason traversible graphs contain a list of source indices for the edges, indicating which edge in the original graph they originate from, and subsequently the index to any associated data.

### Quotient graphs
A quotient graph is a traversible graph where each node represents a set of nodes from the original graph. We use relational indices to create quotient graphs, such as the residue index of atoms.

A quotient of the atom graph to residues will result in each residue being a node, containing all its atoms as constituents. Edges from the original graph will be present only if they connect two of the resulting nodes.

The same function is used for creating basic traversible graphs, subgraphs, and quotient graphs. This may look confusing at first, but comes in handy when combining their functionality.
