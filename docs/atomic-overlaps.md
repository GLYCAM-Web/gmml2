## Atomic overlaps
The main purpose of the carbohydrate and glycoprotein builders is to create 3-dimensional structures of carbohydrates and glycoproteins to be used in molecular dynamics simulations. As such, it is important that
1. Their shapes approximate experimental data, and
2. Interatomic forces are not unreasonably high, as this would destabilize the simulations

The metadata used in gmml2 provides values for bond lengths and relative atomic positions with which residues are internally stable. However, given the relatively high flexibility of glycosidic linkages, setting rotatable bond angles within them to their default values carries a high risk of the molecule folding in a way such that different residues overlap with eachother, resulting in unrealistically high forces between their atoms.

In order to prevent this, we perform an overlap resolution by adjusting the angles of the glycosidic linkages, after all atoms have been created and bonded.

### Overlap measure
We use a modified Lennard-Jones potential to calculate forces between non-bonded atoms, which we use as a measure of the overlap severity. Our main concern is avoiding high repulsive forces, and for this reason we treat the force between atoms whose distance is greater than the sum of their van der Waals radii as zero.

For atoms closer together than the sum of their van der Waals radii, the [van der Waals equation](https://ambermd.org/Questions/vdwequation.pdf) is used. The parameters and van der Waals radii have been sourced from Amber.

We skip calculating this for a few cases:
1. Atoms within the same residue. Metadata is relied on to construct residues with low-energy states internally
2. Atoms directly bonded to eachother. The forces involved in atomic bonds are generally stronger than the forces we are calculating. We abstract this by keeping bond lengths constant
3. Atoms which are both bonded to a third common atom. See above

We always calculate forces between atoms in non-linked residues, and never calculate them for atoms within the same residue. The special case concerns atoms in residues which are linked.
For the sake of simplicity, for any atom belonging to residue A we omit the force calculation with

1. Other atoms directly bonded to residue A, and
2. Any atoms directly bonded to atoms like the above

The data detailing which atoms to skip is stored in a list, and is accessed through the index of the edge linking these residues in the residue graph.

#### Overlap within and between molecules
There are two main ways of calculating overlaps.
One takes a selection and calculates the total overlap of all atoms within said selection. This is an absolute measure, and appropriate to use e.g in the final output.

The other calculates overlap between two selections. This is most useful for comparing the *difference* in overlap severity between multiple possible positions or rotations where the selections are moving relative to eachother, but not relative to themselves.
The return type is a list of doubles, containing the sum of repulsive forces acting on each atom.
### Bounding spheres and intersections
When resolving overlaps, we consider each atom to be a sphere with radius equal to its van der Waals radius. Two atoms are considered overlapping if their spherical representations intersect. Now consider two spheres, A and B, each encompassing the spheres of their own set of atoms. If these spheres do not intersect, then neither may any atom in A intersect any atom in B.

If they do intersect, then only the atoms in B which themselves intersect the sphere A may intersect any atoms in A. Calculating the bounding sphere of each residue lets us use this principle to avoid having to calculate the distance between every single pair of atoms.

In the case of the glycoprotein builder, we extend this by calculating bounding spheres for molecules based on their residues.
The algorithm for finding bounding spheres is based on [this paper](https://www.researchgate.net/publication/242453691_An_Efficient_Bounding_Sphere).

### Linkage shapes and overlap resolution
There are two reasons for modifying the angles of rotatable bonds. First and foremost in order to resolve any overlap, but a secondary reason is to force a particular shape of the glycan. The carbohydrate builder aims to generate shapes using the default angle of the most likely rotamers. The glycoprotein builder allows whichever rotamers give the best fit to be selected, while varying the angles.

When a bond is rotated, the entire molecule can be considered as two sets of atoms, one set on either side of the bond. Atoms belonging to the same set stay in the same position relative to each other when the bond is rotated. Thus, calculating forces for atoms only on opposite sides of the bond is sufficient for comparing the relative effects of different angles on the force.

For each selected rotamer, we generate an even distribution of angles within the allowed range. The increment is given in the parameters. Starting with the preferred angle and working our way out, we calculate rotated positions for the set of atoms on the root side of the bond, then calculate overlap severity between the two sets of atoms.

Whichever angle resulted in the lowest amount of overlap is chosen as the new dihedral angle, and the atom positions updated. In the case where multiple angles resulted in the same lowest overlap, the angle closest to the preferred angle is chosen.

#### Extra searches
Once an angle has been found, the algorithm will make a number of extra searches searches around the best angle, halving the angle increment at each step. This is done in order to improve accuracy without significantly affecting speed. The number of extra searches can be set to zero.

#### Conformers vs permutations
When conformer linkages are concerned, all bonds must use the same rotamer. Thus the order of iteration is switched to iterate over rotamers first, and then over each bond to try different angles using only that rotamer.

The conformer / permutation code handling would be cleaner if permutations were handled as group of conformer linkages with 1 bond each.
