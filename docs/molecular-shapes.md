## Molecular shapes
### Rotatable bonds
Certain atomic bonds can rotate on their axis. Specifically, these are bonds acting as the single point of connection between two parts of a molecule. Given the euclidian positions of the two bonded atoms, `u` and `v`, the axis of rotation is given as the normalized vector `u - v`.

When such a bond is rotated, all the atoms on one side of the bond will spin around the axis, relative to all the atoms on the other side of the bond. For the sake of consistency, we tend to keep the side containing the root node in place, only moving the opposite side. This has the advantage of giving every molecule a fixed point, which helps when comparing variations of the same molecule in a 3d viewer.

### Dihedral angles
In order to assign an angle to a bond, more than two coordinates are needed. For the bond between atoms [b, c], the coordinates of atoms [a, b, c, d] are used, where both [a, b] and [c, d] are bonded. Which atoms these are is defined in metadata.

In some cases, the data is missing hydrogen atoms necessary to calculate certain dihedral angles. When this is occurs, we create the hydrogen atoms, but omit writing them to any output file.
### Glycosidic linkages
Where monosaccarides are concerned, we are mainly interested in rotatable bonds located in the linkages between residues. While rotatable bonds can be found elsewhere, these are located either

* on the fringes of a residue, where the rotation has very little impact on the overall molecular shape
* within a residue, dividing it in two parts with linkages on both sides

Rotatable bonds within monosaccharide residues would be useful to include, but are not supported as of v1.3.

Glycosidic linkages are linkages between two monosaccharides, a monosaccharide and its derivative, or between a monosaccharide and a protein. These linkages contain multiple rotatable bonds (typically 4, but the exceptions are many) with large degrees of rotational freedom.
As such, they are highly flexible, and given the branching character of many carbohydrates, this leads to the same molecule often having a large variety of possible shapes.

Most glycosidic linkages are between two residues only. However, there are Y-shaped linkages in which three residues connect together. The code gets quite convoluted when handling these, and could use being rewritten.
It might also help to redefine the abstraction from “residue linkages” to “groups of rotatable bonds”. Until bonds within monosaccharide residues are handled, the current abstraction is good enough.

### Rotamers
Bonds experience varying amounts of energy at various angles, and as such do not take on all angles equally. Instead of attempting to simulate these forces, gmml2 includes metadata specifying which ranges of angles are considered low in energy. Each such a range is called a rotamer, and a single bond may have multiple rotamers. The range of a rotamer is defined as a default angle, and a freedom of motion either direction. Freedoms of motion are defined either as
1. degrees in one standard deviation
2. the maximum allowed deviation from the default angle, in degrees

As an example, a bond with three rotamers each with a maximum deviation of ±20 degrees may be allowed to take on an angle in the following ranges, but never a value in between
```
[40, 80]
[160, 200]
[280, 320]
```

Glycosidic linkages can be categorized in two types

1. Conformers, in which the rotamers of each bond are dependent on eachother
2. Permutations, in which they are independent

Rotamer ids are currently stored per linkage/bond. We might want to create an intermediate table of bond id, rotamer id pairs, in order to get a better overview and make it easier to modify for special cases.

### Sidechain rotamers
The Dunbrack rotamer library contains data for dihedral angles of up to 4 rotatable bonds per amino acid sidechain. These do not have any range of motion. Each rotamer instead contains a discrete set of angles for each bond in the sidechain. When setting their shape, each angle is applied at once.
The glycoprotein builder makes use of adjusting the shapes of sidechains.
