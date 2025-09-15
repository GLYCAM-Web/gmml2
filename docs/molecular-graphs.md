## Molecular graphs
### Assembly
This section referes to a collection of C++ data structures and functions. It is unrelated to assembly, the programming language. It is also named the same as the macromolecular assembly, which it perhaps should not be.

---

For atomic graphs, the assembly indices data structure provides a useful way of organizing molecular structures into the common atom, residue, molecule arrangement.

It contains lists of indices for
* which residue each atom belongs to
* which molecule each residue belongs to
* which assembly each molecule belongs to

Paired with an atom graph, they contain the data necessary to construct quotient graphs of residues or molecules, with inter-residue connectivity matching that seen in the atom graph.

The assembly code also contains functions for selecting subsets of atoms, residues, or molecules. This is useful when building or using functionality that is not always interested in the entire molecular structure.
### Residue graphs
A useful abstraction of a carbohydrate molecule is the residue graph. This graph (under current assumptions) always takes the form of a tree. That is, it has a single root node and well-defined ordering.

Residue graphs are generated using a text string in the [GLYCAM condensed oligosaccharide notation](https://glycam.org/docs/custombuilders/condensed-notation/index.html). The sequence parser translates from this notation into a tree of residues. Each residue may be either

* an aglycone (which may only be the root)
* a monosaccharide (which may have child nodes connecting to them)
* a derivative/deoxy, which always connect to monosaccharides and may only be leaves.

Residues are defined by their names, and can connect in a number of ways. Monosaccharides connect through positions along their rings defined in positive integers. A monosaccharide always connects to its parent through its 1-position, and may have children connected to any of its 2+ positions. Aglycone, derivatives, and deoxy only have one possible way of connecting. Hence their position is left implicit in the notation. When a position is defined for derivatives and deoxy, it refers to the position on the monosaccharide at which they connect.

#### Heads and tails
The 1-position through which a monosaccharide connects to its parent is called an anomeric linkage. Anomeric linkages are sometimes referred to as *tail*, while the child linkages are referred to as *head*. In other words, the tail is the parent, while the head is the child. An edge would be defined as [tail, head].
### Atomic graphs
While the residue graphs are useful, the main feature of gmml2 is its capability to quickly reproduce experimental data down to the atomic level. To achieve this, we use a combination of residue graphs, parameter files, and overlap resolution algorithms. These are explained in more detail in later sections.
