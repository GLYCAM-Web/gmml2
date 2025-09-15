## Carbohydrates
### Construction
To arrive at an atomic graph of a carbohydrate, a few steps are taken
1. A residue graph is the starting point
2. Data is loaded containing per-residue information of its atoms, their names, elements, charges, connectivity, and relative positions
3. An atomic graph is created, containing entries for each atom of each residue
4. Residue atoms are bonded, and their relative positions adjusted in order to achieve the correct bond lengths
5. Metadata describing likely angles and possible rotations of the bonds of atoms linking residues together is loaded. This is described more in detail in later sections.
6. Atomic overlaps are resolved. This step follows below

### Carbohydrate overlap resolution
Resolving overlaps is part of constructing a carbohydrate. The algorithm here creates what we call a default shape, where we constrain each bond to stay within its most likely rotamer.

Using the described method of calculating overlap severity, the algorithm consists of two steps. First, when the molecule is constructed, the connecting linkage of each residue is adjusted to minimize overlap the moment it is added to the molecule.

When the entire molecule has been constructed, and additional iteration is made over each linkage, starting near the root of the graph.
For each bond, possible angles of the selected rotamer are tried. The angle resulting in the least force is returned. If a force of 0 is found, the function can return early.

The algorithm then proceeds to the next rotatable bond, or the next glycosidic linkage.

When all linkages have been iterated over, the algorithm runs a second time. This has been shown to find shapes with angles closer to those preferred. Continuing the iteration more than two times has not been found to lead to any significant improvement.
