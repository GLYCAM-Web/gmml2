# GMML2 Developer guide
This guide aims to provide an outline of the purpose, capabilities, and interaction of the various parts of the code in [GMML2](https://github.com/GLYCAM-Web/gmml2) from a developer perspective. For learning about specific functions and data structures we refer to the code itself.

### Useful domain knowledge
In molecular structures, atoms are often grouped into residues (e.g monosaccharides or amino acids), and residues grouped into molecules (e.g a carbohydrate or a protein chain).

Molecules can further be assembled into larger structures (e.g multi-chain proteins or glycoproteins), known as macromolecular assemblies. We simply call these assemblies. Most functionality in gmml2 follows this model. The exception is the sequence code, which deals exclusively with residues.

## Key concepts

[Data structures](data-structures.md)

[Graphs](graphs.md)

[Molecular graphs](molecular-graphs.md)

[Molecular shapes](molecular-shapes.md)

[Atomic overlaps](atomic-overlaps.md)

[Carbohydrates](carbohydrates.md)

[Glycoproteins](glycoproteins.md)

[File types](files.md)

### Practical how-to's

[Testing](testing.md)

[Swig wrapping](swig-wrapping.md)

[Releases and versioning](releases-versioning.md)
